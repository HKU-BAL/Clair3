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


/* simple and fast CSV reader:
 * 1. Open CSV file by calling CsvOpen("filename.csv")
 * 2. Read CSV row by calling CsvReadNextRow(csv_handle)
 * 3. Read single CSV line column by calling CsvReadNextCol(returned_row_str, csv_handle)
 */

#ifndef CSV_H_INCLUDED
#define CSV_H_INCLUDED

#ifdef __cplusplus
extern "C" {  /* C++ name mangling */
#endif

/* pointer to private handle structure */
typedef struct CsvHandle_ *CsvHandle;

/**
 * openes csv file
 * @filename: pathname of the file
 * @return: csv handle
 * @notes: you should call CsvClose() to release resources
 */
CsvHandle CsvOpen(const char* filename);
CsvHandle CsvOpen2(const char* filename,
                   char delim,
                   char quote,
                   char escape);

/**
 * closes csv handle, releasing all resources
 * @handle: csv handle
 */
void CsvClose(CsvHandle handle);

/**
 * reads (first / next) line of csv file
 * @handle: csv handle
 */
char* CsvReadNextRow(CsvHandle handle);

/**
 * get column of file
 * @row: csv row (you can use CsvReadNextRow() to parse next line)
 * @context: handle returned by CsvOpen() or CsvOpen2()
 */
const char* CsvReadNextCol(char* row, CsvHandle handle);

#ifdef __cplusplus
};
#endif

#endif
