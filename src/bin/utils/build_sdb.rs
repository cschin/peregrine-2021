// Peregrine Assembler and SHIMMER Genome Assembly Toolkit
// 2019, 2020, 2021- (c) by Jason, Chen-Shan, Chin
//
// This Source Code Form is subject to the terms of the
// Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//
// You should have received a copy of the license along with this
// work. If not, see <http://creativecommons.org/licenses/by-nc-sa/4.0/>.

#![allow(dead_code)]

use flate2::bufread::MultiGzDecoder;
use rayon::prelude::*;
use std::fs::File;
use std::io::prelude::*;
use std::io::{self, BufReader, SeekFrom};

pub struct SeqRec {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
}

enum Fastx {
    FastQ,
    FastA,
}
pub struct FastxReader<R> {
    // struct for reading different file types
    inner: R,
    t: Fastx,
}

impl<R: BufRead> FastxReader<R> {
    pub fn new(mut inner: R, filename: &String) -> Result<Self, io::Error> {
        let t: Fastx;
        {
            // peek the file to decide if it is fasta or fastq
            let r = inner.by_ref();
            let mut buf = Vec::<u8>::new();
            r.take(1).read_to_end(&mut buf)?;
            if buf.len() < 1 {
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    format!("empty file: {}", filename),
                ));
            }
            match buf[0] {
                b'>' => t = Fastx::FastA,
                b'@' => t = Fastx::FastQ,
                _ => t = Fastx::FastA,
            }
        }
        Ok(Self { inner, t })
    }

    pub fn next_rec(&mut self) -> Option<io::Result<SeqRec>> {
        match self.t {
            Fastx::FastA => self.fasta_next_rec(),
            Fastx::FastQ => self.fastq_next_rec(),
        }
    }

    pub fn fasta_next_rec(&mut self) -> Option<io::Result<SeqRec>> {
        // naive partser for fasta format to get the next record

        let mut id_tmp = Vec::<u8>::with_capacity(512);
        let mut seq = Vec::<u8>::with_capacity(1 << 14);

        let res = self.inner.read_until(b'\n', &mut id_tmp);
        if res.is_err() {
            Some(res);
        } else if res.ok() == Some(0) {
            return None;
        }
        let mut r = BufReader::new(&id_tmp[..]);
        let mut id = Vec::<u8>::with_capacity(512);
        let res = r.read_until(b' ', &mut id);
        if res.is_err() {
            Some(res);
        }
        let id = id
            .into_iter()
            .filter(|c| *c != b'\n' && *c != b' ' && *c != b'\r')
            .collect();
        let _x = self.inner.read_until(b'>', &mut seq);
        let seq = seq
            .into_iter()
            .filter(|c| *c != b'\n' && *c != b'>' && *c != b'\r')
            .collect();
        let rec = SeqRec { id: id, seq: seq };

        Some(Ok(rec))
    }

    pub fn fastq_next_rec(&mut self) -> Option<io::Result<SeqRec>> {
        // naive partser for fastq format to get the next record
        // QV strings are ignored

        let mut buf = Vec::<u8>::with_capacity(512);
        let mut id_tmp = Vec::<u8>::with_capacity(512);
        let mut seq = Vec::<u8>::with_capacity(1 << 14);

        let _res = self.inner.read_until(b'\n', &mut id_tmp); //read id
                                                              // fetch the first id up to the first space, strip '\n'
        let mut r = BufReader::new(&id_tmp[..]);
        let mut id = Vec::<u8>::with_capacity(512);
        let _res = r.read_until(b' ', &mut id);
        let id = id
            .into_iter()
            .filter(|c| *c != b'\n' && *c != b' ' && *c != b'\r')
            .collect();
        // get the seq
        let _res = self.inner.read_until(b'\n', &mut seq);
        let seq = seq
            .into_iter()
            .filter(|c| *c != b'\n' && *c != b'\r')
            .collect();
        let rec = SeqRec { id: id, seq: seq };
        // ignore QV
        let _res = self.inner.read_until(b'+', &mut buf);
        let _res = self.inner.read_until(b'\n', &mut buf);
        let _res = self.inner.read_until(b'\n', &mut buf);
        let res = self.inner.read_until(b'@', &mut buf); //get to id line
        if res.is_err() {
            Some(res);
        } else if res.ok() == Some(0) {
            return None;
        }
        Some(Ok(rec))
    }
}

fn get_hpc_flag(seq0: &Vec<u8>) -> Vec<u8> {
    // We use a hybrid approach to handle homopolymer sequence
    // In somce case, it is use to know a base is part of long homopolymer and
    // it is prone to have insertion or deletion errors. We don't compress them
    // when we store the sequences but we mark those bases in case it is useful
    // to ignore them.

    let mut flag = Vec::<u8>::with_capacity(seq0.len());
    let mut i = 0_usize;
    let seq0len = seq0.len();
    while i < seq0len {
        let mut j = i;
        while j < seq0len - 2 && seq0[j] == seq0[j + 1] {
            j += 1;
        }
        if j != i {
            // mask HP > 5 bases
            let mut count = 0_u32;
            while i <= j {
                if count < 4 {
                    flag.push(0b0000);
                } else {
                    flag.push(0b0100);
                }
                count += 1;
                i += 1;
            }
        } else {
            //dimer case
            let mut j = i;
            while j < seq0len - 4 && seq0[j] == seq0[j + 2] && seq0[j + 1] == seq0[j + 3] {
                j += 2;
            }
            if j != i {
                let mut count = 0_u32;
                while i <= j {
                    if count < 4 {
                        flag.push(0b0000);
                        flag.push(0b0000);
                    } else {
                        flag.push(0b1000);
                        flag.push(0b1000);
                    }
                    i += 2;
                    count += 2;
                }
            } else {
                flag.push(0x0000);
                i += 1;
            }
        }
    }
    flag
}

pub fn encode_biseq(s: &Vec<u8>) -> Vec<u8> {
    // we use 4 bits to store the each base, the reversed compliment are stored together
    // the lower 4 bits are for the base from the original orientation
    // the higher 4 bits are for the base from the reversed complementary orientation
    // For each 4 bit field, lower two bits are for the base, the higher two bits are flags.
    // 'A' = 0bxxx00, 'C' = 0bxx01, 'G' = 0bxx10, 'T' = 0bxx11
    // flag: 0b10xx = hp tagged, 0b00xx = non hp tagged
    // 0b1100 (12) : None base

    let fourbit_map_f: [u8; 256] = [
        12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
        12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
        12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 0, 12, 1, 12,
        12, 12, 2, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 3, 12, 12, 12, 12, 12, 12, 12,
        12, 12, 12, 12, 12, 0, 12, 1, 12, 12, 12, 2, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
        12, 3, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
        12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
        12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
        12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
        12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
        12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
        12, 12, 12,
    ];
    let len = s.len();
    let mut out_s = Vec::<u8>::with_capacity(len);
    let mut f_seq = Vec::<u8>::with_capacity(len);
    let mut r_seq = Vec::<u8>::with_capacity(len);
    for p in 0..len {
        let rp = len - 1 - p;
        //let code = ((fourbit_map_f[s[rp] as usize] ^ 0x03) << 4) | fourbit_map_f[s[p] as usize];
        //out_s.push(code);
        f_seq.push(fourbit_map_f[s[p] as usize]);
        r_seq.push(fourbit_map_f[s[rp] as usize] ^ 0b0011);
    }

    let f_flag = get_hpc_flag(&f_seq);
    let r_flag = get_hpc_flag(&r_seq);
    for p in 0..len {
        out_s.push(((r_flag[p] | r_seq[p]) << 4) | (f_flag[p] | f_seq[p]));
    }
    out_s
}

pub fn build(seq_list_file: &String, out_prefix: &String) -> Result<usize, io::Error> {
    // given a list of file in `seq_list_file`, read the sequences to building the sequence
    // database and index

    let seqdb_name = format!("{}.seqdb", out_prefix);
    log::info!("create seq db: {}", seqdb_name);
    let mut out_db_file = File::create(seqdb_name)?;
    let seqidx_name = format!("{}.idx", out_prefix);
    log::info!("create seq index: {}", seqidx_name);
    let mut out_idx_file = File::create(seqidx_name)?;
    let mut start = 0_usize;
    let mut seq_id = 0_u32;

    log::info!("get input files from: {}", seq_list_file);
    let f = File::open(seq_list_file)?;
    let seq_list_buf = BufReader::new(f);

    for fastx_file in seq_list_buf.lines() {
        let input_fn = fastx_file.unwrap();
        log::info!("input file: {}", input_fn);
        let metadata = std::fs::metadata(&input_fn)?;
        if !metadata.is_file() || metadata.len() < (1 << 16) {
            log::info!(
                "input file: {} may not be proper input file (filesize = {}), ignore",
                input_fn,
                metadata.len()
            );
            continue;
        }
        let input_file = File::open(&input_fn)?;
        let mut reader = BufReader::new(input_file);
        let mut is_gzfile = false;
        {
            let r = reader.by_ref();
            let mut buf = Vec::<u8>::new();
            let _ = r.take(2).read_to_end(&mut buf);
            if buf == [0x1F_u8, 0x8B_u8] {
                log::info!("input file: {} detected as gz-compressed file", input_fn);
                is_gzfile = true;
            }
        }

        let _ = reader.seek(SeekFrom::Start(0));
        let mut seqs = Vec::<(u32, Vec<u8>, Vec<u8>)>::new();
        if is_gzfile {
            let fastx_buf = BufReader::new(MultiGzDecoder::new(&mut reader));
            let mut fastx_reader = FastxReader::new(fastx_buf, &input_fn)?;
            while let Some(r) = fastx_reader.next_rec() {
                let r = r.unwrap();
                if r.seq.len() < 500 {
                    //ignore very short reads
                    continue;
                }
                seqs.push((seq_id, r.id, r.seq));
                seq_id += 1;
            }
        } else {
            let mut fastx_reader = FastxReader::new(reader, &input_fn)?;
            while let Some(r) = fastx_reader.next_rec() {
                let r = r.unwrap();
                if r.seq.len() < 500 {
                    //ignore very short reads
                    continue;
                }
                seqs.push((seq_id, r.id, r.seq));
                seq_id += 1;
            }
        }
        let biseq = seqs
            .par_iter()
            .map(|(id, name, s)| (*id, name.clone(), encode_biseq(s)))
            .collect::<Vec<(u32, Vec<u8>, Vec<u8>)>>();

        biseq.iter().for_each(|(id, name, s)| {
            let _ = out_db_file.write(&s);
            let _ = writeln!(
                out_idx_file,
                "{:09} {} {} {}",
                id,
                String::from_utf8_lossy(name),
                s.len(),
                start
            );
            start += s.len();
        });
    }
    log::info!("total number of reads indexed: {}", seq_id);
    log::info!("total number of bases: {}", start);
    log::info!("average read length: {}", start as f32 / seq_id as f32);
    Ok(start)
}
