/* Copyright (c) 2019 Jan Doczy
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "csv.h"

/* Windows specific */
#ifdef _WIN32
#include <Windows.h>
typedef unsigned long long file_off_t;
#else
#include <sys/types.h>
typedef off_t file_off_t;
#endif

/* max allowed buffer */
#define BUFFER_WIDTH_APROX (40 * 1024 * 1024)

#if defined (__aarch64__) || defined (__amd64__) || defined (_M_AMD64)
/* unpack csv newline search */
#define CSV_UNPACK_64_SEARCH
#endif

/* private csv handle:
 * @mem: pointer to memory
 * @pos: position in buffer
 * @size: size of memory chunk
 * @context: context used when processing cols
 * @blockSize: size of mapped block
 * @fileSize: size of opened file
 * @mapSize: ...
 * @auxbuf: auxiliary buffer
 * @auxbufSize: size of aux buffer
 * @auxbufPos: position of aux buffer reader
 * @quotes: number of pending quotes parsed
 * @fh: file handle - descriptor
 * @delim: delimeter - ','
 * @quote: quote '"'
 * @escape: escape char
 */
struct CsvHandle_
{
    void* mem;
    size_t pos;
    size_t size;
    char* context;
    size_t blockSize;
    file_off_t fileSize;
    file_off_t mapSize;
    size_t auxbufSize;
    size_t auxbufPos;
    size_t quotes;
    void* auxbuf;
    
#if defined ( __unix__ )
    int fh;
#elif defined ( _WIN32 )
    HANDLE fh;
    HANDLE fm;
#else
    #error Wrong platform definition
#endif

    char delim;
    char quote;
    char escape;
};

CsvHandle CsvOpen(const char* filename)
{
    /* defaults */
    return CsvOpen2(filename, ',', '"', '\\');
}

/* trivial macro used to get page-aligned buffer size */
#define GET_PAGE_ALIGNED( orig, page ) \
    (((orig) + ((page) - 1)) & ~((page) - 1))

/* thin platform dependent layer so we can use file mapping
 * with winapi and oses following posix specs.
 */
#ifdef __unix__
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

CsvHandle CsvOpen2(const char* filename,
                   char delim,
                   char quote,
                   char escape)
{
    /* alloc zero-initialized mem */
    long pageSize;
    struct stat fs;

    CsvHandle handle = calloc(1, sizeof(struct CsvHandle_));
    if (!handle)
        goto fail;

    /* set chars */
    handle->delim = delim;
    handle->quote = quote;
    handle->escape = escape;

    /* page size */
    pageSize = sysconf(_SC_PAGESIZE);
    if (pageSize < 0)
        goto fail;

    /* align to system page size */
    handle->blockSize = GET_PAGE_ALIGNED(BUFFER_WIDTH_APROX, pageSize);
    
    /* open new fd */
    handle->fh = open(filename, O_RDONLY);
    if (handle->fh < 0)
        goto fail;

    /* get real file size */
    if (fstat(handle->fh, &fs))
    {
       close(handle->fh);
       goto fail;
    }
    
    handle->fileSize = fs.st_size;
    return handle;
    
  fail:
    free(handle);
    return NULL;
}

static void* MapMem(CsvHandle handle)
{
    handle->mem = mmap(0, handle->blockSize,
                       PROT_READ | PROT_WRITE,
                       MAP_PRIVATE,
                       handle->fh, handle->mapSize);
    return handle->mem;
}

static void UnmapMem(CsvHandle handle)
{
    if (handle->mem)
        munmap(handle->mem, handle->blockSize);
}

void CsvClose(CsvHandle handle)
{
    if (!handle)
        return;

    UnmapMem(handle);

    close(handle->fh);
    free(handle->auxbuf);
    free(handle);
}

#else

/* extra Windows specific implementations
 */
CsvHandle CsvOpen2(const char* filename,
                   char delim,
                   char quote,
                   char escape)
{
    LARGE_INTEGER fsize;
    SYSTEM_INFO info;
    size_t pageSize = 0;
    CsvHandle handle = calloc(1, sizeof(struct CsvHandle_));
    if (!handle)
        return NULL;

    handle->delim = delim;
    handle->quote = quote;
    handle->escape = escape;

    GetSystemInfo(&info);
    handle->blockSize = GET_PAGE_ALIGNED(BUFFER_WIDTH_APROX, info.dwPageSize);
    handle->fh = CreateFile(filename, 
                            GENERIC_READ, 
                            FILE_SHARE_READ, 
                            NULL, 
                            OPEN_EXISTING, 
                            FILE_ATTRIBUTE_NORMAL, 
                            NULL);

    if (handle->fh == INVALID_HANDLE_VALUE)
        goto fail;

    if (GetFileSizeEx(handle->fh, &fsize) == FALSE)
        goto fail;

    handle->fileSize = fsize.QuadPart;
    if (!handle->fileSize)
        goto fail;

    handle->fm = CreateFileMapping(handle->fh, NULL, PAGE_WRITECOPY, 0, 0, NULL);
    if (handle->fm == NULL)
        goto fail;

    return handle;

fail:
    if (handle->fh != INVALID_HANDLE_VALUE)
        CloseHandle(handle->fh);

    free(handle);
    return NULL;
}

static void* MapMem(CsvHandle handle)
{
    size_t size = handle->blockSize;
    if (handle->mapSize + size > handle->fileSize)
        size = 0;  /* last chunk, extend to file mapping max */

    handle->mem = MapViewOfFileEx(handle->fm, 
                                  FILE_MAP_COPY,
                                  (DWORD)(handle->mapSize >> 32),
                                  (DWORD)(handle->mapSize & 0xFFFFFFFF),
                                  size,
                                  NULL);
    return handle->mem;
}

static void UnmapMem(CsvHandle handle)
{
    if (handle->mem)
        UnmapViewOfFileEx(handle->mem, 0);
}

void CsvClose(CsvHandle handle)
{
    if (!handle)
        return;

    UnmapMem(handle);

    CloseHandle(handle->fm);
    CloseHandle(handle->fh);
    free(handle->auxbuf);
    free(handle);
}

#endif

static int CsvEnsureMapped(CsvHandle handle)
{
    file_off_t newSize;
    
    /* do not need to map */
    if (handle->pos < handle->size)
        return 0;

    UnmapMem(handle);  

    handle->mem = NULL;
    if (handle->mapSize >= handle->fileSize)
        return -EINVAL;

    newSize = handle->mapSize + handle->blockSize;
    if (MapMem(handle))
    {
        handle->pos = 0;
        handle->mapSize = newSize;

        /* read only up to filesize:
         * 1. mapped block size is < then filesize: (use blocksize)
         * 2. mapped block size is > then filesize: (use remaining filesize) */
        handle->size = handle->blockSize;
        if (handle->mapSize > handle->fileSize)
            handle->size = (size_t)(handle->fileSize % handle->blockSize);
        
        return 0;
    }
    
    return -ENOMEM;
}

static char* CsvChunkToAuxBuf(CsvHandle handle, char* p, size_t size)
{
    size_t newSize = handle->auxbufPos + size + 1;
    if (handle->auxbufSize < newSize)
    {
        void* mem = realloc(handle->auxbuf, newSize);
        if (!mem)
            return NULL;

        handle->auxbuf = mem;
        handle->auxbufSize = newSize;
    }

    memcpy((char*)handle->auxbuf + handle->auxbufPos, p, size);
    handle->auxbufPos += size;
    
    *(char*)((char*)handle->auxbuf + handle->auxbufPos) = '\0';
    return handle->auxbuf;
}

static void CsvTerminateLine(char* p, size_t size)
{
    /* we do support standard POSIX LF sequence
     * and Windows CR LF sequence.
     * old non POSIX Mac OS CR is not supported.
     */
    char* res = p;
    if (size >= 2 && p[-1] == '\r')
        --res;

    *res = 0;
}

#define CSV_QUOTE_BR(c, n) \
    do \
        if (c##n == quote)                              \
            handle->quotes++;                           \
        else if (c##n == '\n' && !(handle->quotes & 1)) \
            return p + n;                               \
    while (0)


static char* CsvSearchLf(char* p, size_t size, CsvHandle handle)
{
    /* TODO: this can be greatly optimized by
     * using modern SIMD instructions, but for now
     * we only fetch 8Bytes "at once"
     */
    char* end = p + size;
    char quote = handle->quote;

#ifdef CSV_UNPACK_64_SEARCH
    uint64_t* pd = (uint64_t*)p;
    uint64_t* pde = pd + (size / sizeof(uint64_t));

    for (; pd < pde; pd++)
    {
        /* unpack 64bits to 8x8bits */
        char c0, c1, c2, c3, c4, c5, c6, c7;
        p = (char*)pd;
        c0 = p[0];
        c1 = p[1];
        c2 = p[2];
        c3 = p[3];
        c4 = p[4];
        c5 = p[5];
        c6 = p[6];
        c7 = p[7];

        CSV_QUOTE_BR(c, 0);
        CSV_QUOTE_BR(c, 1);
        CSV_QUOTE_BR(c, 2);
        CSV_QUOTE_BR(c, 3);
        CSV_QUOTE_BR(c, 4);
        CSV_QUOTE_BR(c, 5);
        CSV_QUOTE_BR(c, 6);
        CSV_QUOTE_BR(c, 7);
    }
    p = (char*)pde;
#endif
    
    for (; p < end; p++)
    {
        char c0 = *p;
        CSV_QUOTE_BR(c, 0);
    }

    return NULL;
}

char* CsvReadNextRow(CsvHandle handle)
{
    size_t size;
    char* p = NULL;
    char* found = NULL;

    do
    {
        int err = CsvEnsureMapped(handle);
        handle->context = NULL;
        
        if (err == -EINVAL)
        {
            /* if this is n-th iteration
             * return auxbuf (remaining bytes of the file) */
            if (p == NULL)
                break;

            return handle->auxbuf;
        }
        else if (err == -ENOMEM)
        {
            break;
        }
        
        size = handle->size - handle->pos;
        if (!size)
            break;

        /* search this chunk for NL */
        p = (char*)handle->mem + handle->pos;
        found = CsvSearchLf(p, size, handle);

        if (found)
        {
            /* prepare position for next iteration */
            size = (size_t)(found - p) + 1;
            handle->pos += size;
            handle->quotes = 0;
            
            if (handle->auxbufPos)
            {
                if (!CsvChunkToAuxBuf(handle, p, size))
                    break;
                
                p = handle->auxbuf;
                size = handle->auxbufPos;
            }

            /* reset auxbuf position */
            handle->auxbufPos = 0;

            /* terminate line */
            CsvTerminateLine(p + size - 1, size);
            return p;
        }
        else
        {
            /* reset on next iteration */
            handle->pos = handle->size;
        }

        /* correctly process boundries, storing
         * remaning bytes in aux buffer */
        if (!CsvChunkToAuxBuf(handle, p, size))
            break;

    } while (!found);

    return NULL;
}

const char* CsvReadNextCol(char* row, CsvHandle handle)
{
    /* return properly escaped CSV col
     * RFC: [https://tools.ietf.org/html/rfc4180]
     */
    char* p = handle->context ? handle->context : row;
    char* d = p; /* destination */
    char* b = p; /* begin */
    int quoted = 0; /* idicates quoted string */

    quoted = *p == handle->quote;
    if (quoted)
        p++;

    for (; *p; p++, d++)
    {
        /* double quote is present if (1) */
        int dq = 0;
        
        /* skip escape */
        if (*p == handle->escape && p[1])
            p++;

        /* skip double-quote */
        if (*p == handle->quote && p[1] == handle->quote)
        {
            dq = 1;
            p++;
        }

        /* check if we should end */
        if (quoted && !dq)
        {
            if (*p == handle->quote)
                break;
        }
        else if (*p == handle->delim)
        {
            break;
        }

        /* copy if required */
        if (d != p)
            *d = *p;
    }
    
    if (!*p)
    {
        /* nothing to do */
        if (p == b)
            return NULL;

        handle->context = p;
    }
    else
    {
        /* end reached, skip */
        *d = '\0';
        if (quoted)
        {
            for (p++; *p; p++)
                if (*p == handle->delim)
                    break;

            if (*p)
                p++;
            
            handle->context = p;
        }
        else
        {
            handle->context = p + 1;
        }
    }
    return b;
}
