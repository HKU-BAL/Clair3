use anyhow::{Result, bail};
use rust_htslib::{
    bam::{IndexedReader, Read as BamRead, Record},
    htslib,
};

pub(crate) struct Pileup<'a> {
    pub position: i32,
    pub alignments: &'a [htslib::bam_pileup1_t],
}

pub(crate) struct Pileups<'a> {
    reader: &'a mut IndexedReader,
    inner: htslib::bam_plp_t,
}

extern "C" fn next_record(data: *mut libc::c_void, record: *mut htslib::bam1_t) -> i32 {
    let reader = unsafe { &mut *(data as *mut IndexedReader) };
    let mut record_wrapper = Record::new();
    match reader.read(&mut record_wrapper) {
        Some(Ok(())) => {
            // SAFETY: HTSlib owns the transferred allocation; Rust drops only the header.
            unsafe { *record = record_wrapper.inner };
            record_wrapper.inner.data = std::ptr::null_mut();
            0
        }
        _ => -1,
    }
}

impl<'a> Pileups<'a> {
    pub fn new(reader: &'a mut IndexedReader) -> Result<Self> {
        let inner = unsafe {
            htslib::bam_plp_init(Some(next_record), reader as *mut _ as *mut libc::c_void)
        };
        if inner.is_null() {
            bail!("failed to initialize BAM pileup");
        }
        Ok(Self { reader, inner })
    }

    pub fn set_max_depth(&mut self, depth: i32) {
        unsafe { htslib::bam_plp_set_maxcnt(self.inner, depth) };
    }
}

impl Iterator for Pileups<'_> {
    type Item = Result<Pileup<'static>>;

    fn next(&mut self) -> Option<Self::Item> {
        let (mut tid, mut position, mut depth) = (0, 0, 0);
        let alignments =
            unsafe { htslib::bam_plp_auto(self.inner, &mut tid, &mut position, &mut depth) };
        if alignments.is_null() {
            return if depth == -1 {
                Some(Err(anyhow::anyhow!("failed to read BAM pileup")))
            } else {
                None
            };
        }
        // SAFETY: HTSlib owns this alignment array until the next pileup call.
        let alignments = unsafe { std::slice::from_raw_parts(alignments, depth as usize) };
        Some(Ok(Pileup {
            position,
            alignments,
        }))
    }
}

impl Drop for Pileups<'_> {
    fn drop(&mut self) {
        let _ = self.reader;
        unsafe { htslib::bam_plp_destroy(self.inner) };
    }
}

#[derive(Clone, Copy)]
pub(crate) struct RawRecord<'a> {
    inner: &'a htslib::bam1_t,
}

impl<'a> RawRecord<'a> {
    pub fn new(inner: *const htslib::bam1_t) -> Self {
        Self {
            // SAFETY: pileup alignments own their records for the duration of this view.
            inner: unsafe { &*inner },
        }
    }

    pub fn flags(&self) -> u16 {
        self.inner.core.flag
    }

    pub fn mapq(&self) -> u8 {
        self.inner.core.qual
    }

    pub fn pos(&self) -> i64 {
        self.inner.core.pos
    }

    pub fn seq_len(&self) -> usize {
        self.inner.core.l_qseq as usize
    }

    fn sequence(&self) -> &[u8] {
        let offset = self.inner.core.l_qname as usize + self.inner.core.n_cigar as usize * 4;
        let length = (self.seq_len() + 1) / 2;
        // SAFETY: BAM stores packed sequence after qname and CIGAR, within l_data.
        unsafe { std::slice::from_raw_parts(self.inner.data.add(offset) as *const u8, length) }
    }

    pub fn encoded_base(&self, index: usize) -> u8 {
        let encoded = self.sequence()[index / 2];
        (encoded >> ((!index & 1) << 2)) & 0x0f
    }

    pub fn cigar(&self) -> &[u32] {
        let count = self.inner.core.n_cigar as usize;
        if count == 0 {
            return &[];
        }
        // SAFETY: BAM stores n_cigar 32-bit operations immediately after qname.
        unsafe {
            std::slice::from_raw_parts(
                self.inner.data.add(self.inner.core.l_qname as usize) as *const u32,
                count,
            )
        }
    }
}
