// Peregrine Assembler and SHIMMER Genome Assembly Toolkit 
// 2019, 2020, 2021- (c) by Jason, Chen-Shan, Chin
//
// This Source Code Form is subject to the terms of the 
// Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//
// You should have received a copy of the license along with this
// work. If not, see <http://creativecommons.org/licenses/by-nc-sa/4.0/>.

#![allow(dead_code)]

//
// generate error corrected reads from the overlaps
//

use super::shmmrutils::*;
use super::Parameters;
use glob::glob;
use memmap::{Mmap, MmapOptions};
use rustc_hash::FxHashMap;
use rustc_hash::FxHashSet;
use std::fs::File;
use std::io::{self, BufRead, Write};
use std::path::Path;
use std::str::from_utf8;
use std::thread;

use std::io::prelude::*;
use threadpool::ThreadPool;

#[derive(Clone, Copy)]
struct ReadLocation {
    start: usize,
    len: usize,
}

#[derive(Debug, Copy, Clone)]
struct Overlap {
    rid0: u32,
    rid1: u32,
    strand1: u8,
    len0: u32,
    len1: u32,
    d_left: i32,
    d_right: i32,
    bgn0: u32,
    end0: u32,
    bgn1: u32,
    end1: u32,
    dist: u32,
    idt: f32,
    dist_c: u32,
    max_dist_c: u32,
    idt_c: f32,
    flag: u8,
}
// flag bit field
// 0x01: the rid0 is chimer
// 0x02: the rir0 and rid1 are compatitable pair
// 0x04: the rid1 is the best right pair
// 0x08: the rid1 is the best left pair
// 0x10: the rid1 is a chimer
// 0x20: the rid1 is contained
// 0x40: the rid0 is contained

fn build_overlap(v: Vec<&str>) -> Overlap {
    Overlap {
        rid0: v[1].parse().unwrap(),
        rid1: v[2].parse().unwrap(),
        strand1: v[3].parse().unwrap(),
        len0: v[4].parse().unwrap(),
        len1: v[5].parse().unwrap(),
        d_left: v[6].parse().unwrap(),
        d_right: v[7].parse().unwrap(),
        bgn0: v[8].parse().unwrap(),
        end0: v[9].parse().unwrap(),
        bgn1: v[10].parse().unwrap(),
        end1: v[11].parse().unwrap(),
        dist: v[12].parse().unwrap(),
        idt: v[13].parse().unwrap(),
        dist_c: v[14].parse().unwrap(),
        max_dist_c: v[15].parse().unwrap(),
        idt_c: v[16].parse().unwrap(),
        flag: v[17].parse().unwrap(),
    }
}

fn _format_overlap(o: Overlap) -> String {
    format!(
        "{} {} {} {} {} {} {} {} {} {} {} {} {:.2} {} {} {:.2} {}",
        o.rid0,
        o.rid1,
        o.strand1,
        o.len0,
        o.len1,
        o.d_left,
        o.d_right,
        o.bgn0,
        o.end0,
        o.bgn1,
        o.end1,
        o.dist,
        o.idt,
        o.dist_c,
        o.max_dist_c,
        o.idt_c,
        o.flag
    )
}

type OverlapMap = FxHashMap<u32, Vec<Overlap>>;
type _DeltaMap = FxHashMap<(u32, u32, u8, i32), Vec<u32>>;

fn build_read_ovlp_data<P>(filename: P) -> Result<OverlapMap, std::io::Error>
where
    P: AsRef<Path>,
{
    let mut rid2ovlp = OverlapMap::default();
    let mut buffer = String::new();

    let mut file = File::open(filename)?;
    file.read_to_string(&mut buffer)?;
    for line in buffer.split("\n") {
        let mut v: Vec<&str> = Vec::<&str>::with_capacity(24); // we need pre-allocate some space for performance
        line.split(' ').for_each(|c| v.push(c));
        match v[0] {
            "O" => {
                let ovlp = build_overlap(v);
                //d_left = ovlp.d_left;
                rid2ovlp
                    .entry(ovlp.rid0)
                    .or_insert_with(|| vec![])
                    .push(ovlp);
            }
            _ => (),
        }
    }
    Ok(rid2ovlp)
}

fn get_consensue_seq(
    seq0: Vec<u8>,
    support_seq: Vec<([usize; 4], Vec<u8>)>,
    tol: f64,
    min_ec_cov: Option<u16>,
) -> Option<Vec<u8>> {
    let min_ec_cov = min_ec_cov.unwrap_or(1);
    // 4bit encoded seq consensus
    let mut cov = vec![0 as u16; seq0.len()];
    let mut deltas = FxHashMap::<u32, u16>::default();
    for (offsets, seq1) in support_seq {
        let [b0, e0, b1, e1] = offsets;
        let seq0_trim = seq0[b0..e0].to_vec();
        let seq1_trim = seq1[b1..e1].to_vec();
        if let Some(ovlpmatch) = match_reads(&seq0_trim, &seq1_trim, true, tol, 1200, 32) {
            for p in ovlpmatch.bgn0..ovlpmatch.end0 {
                cov[(p as usize) + b0] += 1;
            }
            let mut dpts = ovlpmatch.deltas.unwrap();
            dpts.reverse();
            let mut d = 0_u8;
            let mut px = 0_u32;

            for dpt in dpts {
                let c = if dpt.dk > 0 {
                    0
                } else {
                    seq1_trim[dpt.y as usize - 1]
                };
                let cx = dpt.x + (b0 as u32) - 1;
                if cx != px {
                    d = 0;
                } else {
                    d += 1;
                }

                let key = cx << 12 | (d as u32) << 4 | (c as u32); //this limit the length of the read to be corrected to 2^20
                let counter = deltas.entry(key).or_insert(0);
                *counter += 1;

                px = cx;
            }
        }
    }

    let mut delta_best = FxHashMap::<u32, (u16, u8)>::default();
    for k in deltas.keys() {
        let v = deltas.get(k).unwrap();
        let p = k >> 12;
        if *v < cov[p as usize] >> 1 || *v < 3 {
            continue;
        }
        let key = *k >> 4;
        let counter = delta_best.entry(key).or_insert((*v, (*k & 0xF) as u8));
        if *v > counter.0 {
            delta_best.insert(key, (*v, (*k & 0xF) as u8));
        }
    }

    let mut keys = delta_best.keys().collect::<Vec<&u32>>();
    keys.sort();
    let mut max_delta = FxHashMap::<u32, u8>::default();
    for k in keys.iter() {
        //let v = delta_best.get(k).unwrap();
        let x = *k >> 8;
        max_delta.insert(x, (*k & 0xFF) as u8);
    }

    let mut consensus_seq = Vec::<u8>::new();
    let mut consensus_seq_cov = Vec::<u16>::new();
    for p in 0..seq0.len() {
        if !max_delta.contains_key(&(p as u32)) {
            consensus_seq.push(seq0[p]);
            consensus_seq_cov.push(cov[p]);
            continue;
        }
        let max_d = *max_delta.get(&(p as u32)).unwrap();
        for d in 0..max_d + 1 {
            let k = ((p as u32) << 8) | (d as u32);
            if let Some(v) = delta_best.get(&k) {
                if v.1 != 0 {
                    if d == 0 {
                        consensus_seq.push(seq0[p]);
                        consensus_seq_cov.push(cov[p]);
                    }
                    consensus_seq.push(v.1);
                    consensus_seq_cov.push(cov[p]);
                    // println!("insert {} {}", p, b);
                }
            }
        }
    }
    assert!(consensus_seq_cov.len() == consensus_seq.len());

    let mut bgn = Option::<usize>::None;
    for i in 0..consensus_seq_cov.len() {
        if consensus_seq_cov[i] >= min_ec_cov {
            // need two other supporting reads for EC
            bgn = Some(i);
            break;
        }
    }

    //consensus_seq_cov.iter().for_each(|c| {println!("{}", c)});

    let mut end = Option::<usize>::None;
    let len = consensus_seq_cov.len() - 1;
    for i in 0..consensus_seq_cov.len() {
        if consensus_seq_cov[len - i] >= min_ec_cov {
            end = Some(len - i);
            break;
        }
    }

    if bgn.is_none() || end.is_none() {
        None
    } else {
        let b = bgn.unwrap();
        let e = end.unwrap();
        if min_ec_cov >= 3 {
            for i in b..e {
                if consensus_seq_cov[i] < min_ec_cov {
                    return None;
                }
            }
        }
        if b > e || e < b + 1000 {
            None
        } else {
            Some(consensus_seq[b..e].to_vec())
        }
    }
}

fn get_corrected_seq(
    r: u32,
    ovlp: &Vec<Overlap>,
    read_index: &Vec<ReadLocation>,
    readsdb: &Mmap,
    tol: f64,
    min_ec_cov: Option<u16>,
) -> Option<Vec<u8>> {
    let rid0 = r;
    let rloc0 = read_index[rid0 as usize];
    let s0 = rloc0.start;
    let len0 = rloc0.len;
    let e0 = s0 + len0;
    let mut seq0 = Vec::<u8>::with_capacity(20480);

    let basemap = [0b0001, 0b0010, 0b0100, 0b1000];
    for c in &readsdb[s0..e0] {
        seq0.push(basemap[(c & 0b0011) as usize]);
    }

    // println!("seq {} {}", r, seq2string(&seq0));
    let mut support_seqs = Vec::<([usize; 4], Vec<u8>)>::new();
    for vv in ovlp.iter() {
        if vv.flag & 0b0010 == 0 {
            continue;
        }
        let rid1 = vv.rid1;
        let rloc1 = read_index[rid1 as usize];
        let strand: u8 = vv.strand1;
        let s1 = rloc1.start;
        let len1 = rloc1.len;
        let e1 = s1 + len1;
        let mut seq1 = Vec::<u8>::with_capacity(20480);

        // we need to map 2bit to 4bit so we can use 0x0000 for special cases
        if strand == 0 {
            for c in &readsdb[s1..e1] {
                // seq1.push(basemap[(c & 0x0F) as usize]);
                seq1.push(basemap[(c & 0b0011) as usize]);
            }
        } else {
            for c in &readsdb[s1..e1] {
                //seq1.push(basemap[((c >> 4) & 0x0F) as usize]);
                seq1.push(basemap[((c >> 4) & 0b0011) as usize]);
            }
        }
        let b0 = vv.bgn0 as usize;
        let e0 = vv.end0 as usize;
        let b1 = vv.bgn1 as usize;
        let e1 = vv.end1 as usize;
        if b1 > e1 {
            continue;
        }
        support_seqs.push(([b0, e0, b1, e1], seq1));
    }
    if let Some(c_seq) = get_consensue_seq(seq0, support_seqs, tol, min_ec_cov) {
        let mut out_seq = Vec::<u8>::new();
        for c in c_seq {
            match c {
                0b0001 => out_seq.push(b'A'),
                0b0010 => out_seq.push(b'C'),
                0b0100 => out_seq.push(b'G'),
                0b1000 => out_seq.push(b'T'),
                _ => out_seq.push(b'N'),
            }
        }
        Some(out_seq)
    } else {
        None
    }
}

pub fn ovlp_ec(
    seqdb_file: &String,
    index_file: &String,
    prefix: &String,
    out_prefix: &String,
    parameters: &Parameters,
) -> Result<(), io::Error> {
    let mut read_index = Vec::<ReadLocation>::new();
    let mut read_name = Vec::<String>::new();

    log::info!("read seq idx: {}", seqdb_file);
    let lines = read_lines(index_file)?;
    for line in lines {
        if let Ok(rec) = line {
            // the record line looks like 000000023 m64062_190803_042216/144/ccs 20359 467415
            let v: Vec<&str> = rec.split_whitespace().collect();
            // let rid: u32 = v[0].parse().unwrap();
            let start: usize = v[3].parse().unwrap();
            let len: usize = v[2].parse().unwrap();
            read_index.push(ReadLocation {
                start: start,
                len: len,
            });
            read_name.push(v[1].to_string());
            //println!("{} {} {}", rid, start, len);
        }
    }

    let read_index = read_index;
    log::info!("read seq db: {}", seqdb_file);
    let file = File::open(seqdb_file)?;
    let readsdb = unsafe { MmapOptions::new().map(&file)? };
    //println!("{:?}", args);
    let _ = log::info!("prefix: {}", prefix);
    let _ = log::info!("out_prefix: {}", out_prefix);
    let infile_pattern = [prefix.clone(), "*".to_string()].concat();

    let mut children = Vec::new();
    //let mut chunk: u8 = 0;
    for entry in glob(&infile_pattern).expect("Failed to read glob pattern") {
        match entry {
            Ok(path) => {
                // println!("{:?}", path.display());
                let child = thread::spawn(move || {
                    let rid2ovlp = build_read_ovlp_data(path);
                    rid2ovlp
                });
                children.push(child);
            }
            Err(e) => log::error!("{:?}", e),
        }
        //chunk += 1;
    }

    let mut rid2ovlp_all = OverlapMap::default();
    for child in children {
        let rid2ovlp_p = child.join().expect("oops! the child thread panicked")?;
        rid2ovlp_all.extend(rid2ovlp_p);
    }

    let contained = FxHashSet::<u32>::default();

    let readsdb = std::sync::Arc::new(readsdb);
    let read_index = std::sync::Arc::new(read_index);
    let contained = std::sync::Arc::new(contained);
    let rid2ovlp_all = std::sync::Arc::new(rid2ovlp_all);
    let read_name = std::sync::Arc::new(read_name);
    let pool = ThreadPool::new(parameters.nthreads as usize);
    for chunk in 0..parameters.nchunks {
        let readsdb = readsdb.clone();
        let read_index = read_index.clone();
        let contained = contained.clone();
        let rid2ovlp_all = rid2ovlp_all.clone();
        let read_name = read_name.clone();
        let out_prefix = out_prefix.clone();
        let nchunks = parameters.nchunks;
        let tol = parameters.tol;
        let min_ec_cov = parameters.min_ec_cov;
        pool.execute(move || {
            let out_file_name = format!("{}_{:02}.fa", out_prefix, chunk + 1);
            let mut f = File::create(out_file_name).unwrap();
            for r in rid2ovlp_all.keys() {
                if (*r % nchunks) != chunk {
                    continue;
                }
                if contained.contains(r) {
                    continue;
                }
                let ovlp = rid2ovlp_all.get(r).unwrap();
                // let mut deltas = HashMap::<u32, u8>::new();

                let rid0 = *r;
                if let Some(corrected_seq) =
                    get_corrected_seq(rid0, &ovlp, &read_index, &readsdb, tol, Some(min_ec_cov))
                {
                    let _ = writeln!(
                        f,
                        ">{} {}\n{}",
                        read_name[*r as usize],
                        r,
                        from_utf8(&corrected_seq).unwrap()
                    );
                }
            }
        });
    }
    pool.join();
    Ok(())
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::with_capacity(10000000, file).lines())
}
