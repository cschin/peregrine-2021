// Peregrine Assembler and SHIMMER Genome Assembly Toolkit 
// 2019, 2020, 2021- (c) by Jason, Chen-Shan, Chin
//
// This Source Code Form is subject to the terms of the 
// Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//
// You should have received a copy of the license along with this
// work. If not, see <http://creativecommons.org/licenses/by-nc-sa/4.0/>.


#![allow(dead_code)]

use super::shmmrutils::sequence_to_shmmrs;
use super::shmmrutils::{get_2bit_fragment, ReadLocation};
use super::Parameters;
use byteorder::{LittleEndian, WriteBytesExt};
use memmap::{Mmap, MmapOptions};
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write};
use std::mem::size_of;
use std::path::Path;
use threadpool::ThreadPool;

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn index_chunk(
    chunk: u32,
    total_chunk: u32,
    readsdb: &Mmap,
    read_index: &Vec<ReadLocation>,
    prefix: &String,
    wsize: u32,
    ksize: u32,
    rfactor: u32,
) -> Result<(), io::Error> {
    // Create index for a chunk from the read database
    
    let filename = format!("{}-{:02}-of-{:02}.dat", prefix, chunk, total_chunk);
    let mut out_f = BufWriter::new(File::create(filename).unwrap());

    let mut wrt = Vec::<u8>::with_capacity(1 << 16);
    for seq_id in 0..read_index.len() {
        if seq_id % (total_chunk as usize) != (chunk % total_chunk) as usize {
            continue;
        }
        let len = read_index[seq_id].len as u32;
        log::debug!("build_idx: len: {} {}", seq_id, len);
        let seq = get_2bit_fragment(seq_id as u32, 0, 0, len, &readsdb, read_index);
        let shmmrs = sequence_to_shmmrs(seq_id as u32, &seq, wsize, ksize, rfactor);
        for m in shmmrs {
            wrt.write_u64::<LittleEndian>(m.x)?;
            wrt.write_u64::<LittleEndian>(m.y)?;
        }
    }
    let us = size_of::<usize>();
    assert!(us == 8 as usize); //make sure the usize is a 64bit int.
    out_f.write_u64::<LittleEndian>((wrt.len() >> 4) as u64)?;

    // Not sure if it a bug in BufferWriter, it can not write more than 2Gb at once
    // we chop up the data and write in small chunks
    let c = 24_usize;
    for i in 0..=wrt.len() >> c {
        if ((i + 1) << c) < wrt.len() {
            let s = i << c;
            let e = (i + 1) << c;
            out_f.write(&wrt[s..e])?;
        } else {
            let s = i << c;
            let e = wrt.len();
            out_f.write(&wrt[s..e])?;
        }
    }
    out_f.flush()?;
    Ok(())
}

pub fn build(
    seqdb_file: &String,
    index_file: &String,
    out_prefix: &String,
    parameters: &Parameters,
) -> () {
    // Using thread pool to build the SHIMMER index in paralle
    
    let mut read_index = Vec::<ReadLocation>::new();

    if let Ok(lines) = read_lines(index_file) {
        for line in lines {
            if let Ok(rec) = line {
                let v: Vec<&str> = rec.split_whitespace().collect();
                let start: usize = v[3].parse().unwrap();
                let len: usize = v[2].parse().unwrap();
                read_index.push(ReadLocation {
                    start: start,
                    len: len,
                });
            }
        }
    }
    let mmap_seqdb = File::open(seqdb_file).unwrap();
    let mmap_seqdb = unsafe { MmapOptions::new().map(&mmap_seqdb).unwrap() };

    let read_index = std::sync::Arc::new(read_index);
    let mmap_seqdb = std::sync::Arc::new(mmap_seqdb);

    let pool = ThreadPool::new(parameters.nthreads as usize);

    let _nchunks = parameters.nchunks;
    for i in 0.._nchunks {
        let mmap_seqdb = mmap_seqdb.clone();
        let read_index = read_index.clone();
        let out_prefix = out_prefix.clone();
        let parameters = (*parameters).clone();
        pool.execute(move || {
            let r = index_chunk(
                i + 1,
                _nchunks,
                &mmap_seqdb,
                &read_index,
                &out_prefix,
                parameters.w,
                parameters.k,
                parameters.r,
            );
            match r {
                Err(error) => panic!("build index fail: {}", error),
                Ok(()) => (),
            };
        });
    }
    pool.join();
}
