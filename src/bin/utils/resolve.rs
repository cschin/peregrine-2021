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
// for resolve contigs that are highly similar to each others which are most
// likily be homologuous pairs in a diploid genome
//
use super::build_sdb::FastxReader;
use super::layout::log_asm_summary;
use super::shmmrutils::{sequence_to_shmmrs, MM128};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

use rustc_hash::FxHashMap;

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn base2twobit(s: &Vec<u8>) -> Vec<u8> {
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
    for p in 0..len {
        out_s.push(fourbit_map_f[s[p] as usize]);
    }
    out_s
}

pub fn resolve_ht(
    fasta_file: &String,
    output_prefix: &String,
    w: u32,
    k: u32,
    r: u32,
) -> Result<(), io::Error> {
    log::info!("resolve:fasta_file: {}", fasta_file);
    log::info!("resolve:output_prefix: {}", output_prefix);
    log::info!("resolve:parameters: w:{}, k:{}, r:{}", w, k, r,);

    let input_file = File::open(&fasta_file)?;
    let reader = BufReader::new(input_file);
    let mut fastx_reader = FastxReader::new(reader, &fasta_file)?;

    //let mut seqdb = FxHashMap::<Vec<u8>, Vec<u8>>::default();
    let mut shmmr_db = FxHashMap::<String, (u32, Vec<MM128>)>::default();
    let mut shmmr_map = FxHashMap::<u64, Vec<u64>>::default();
    let mut id2name = FxHashMap::<u32, String>::default();
    let mut rid = 0;
    let mut seq_db = FxHashMap::<String, Vec<u8>>::default();
    while let Some(rec) = fastx_reader.next_rec() {
        let rec = rec.unwrap();
        //seqdb.insert(r.id, r.seq);
        //println!("N {}", String::from_utf8_lossy(&r.id));
        let rec_2bitseq = base2twobit(&rec.seq);
        let shmmers = sequence_to_shmmrs(rid, &rec_2bitseq, w, k, r);
        for mm in shmmers.iter() {
            let hash = mm.x >> 8;
            shmmr_map.entry(hash).or_insert_with(|| vec![]).push(mm.y);
        }
        let n = String::from_utf8_lossy(&rec.id).into_owned();
        id2name.insert(rid, n.clone());
        shmmr_db.insert(n.clone(), (rid, shmmers));
        seq_db.insert(n, rec.seq);
        rid += 1;
    }

    let mut matches = FxHashMap::<(u32, u32), Vec<(u32, u32)>>::default();

    for (_, (rid, shmmrs)) in shmmr_db.iter() {
        for mer0 in shmmrs {
            let y0 = mer0.y;
            let pos0 = ((y0 & 0xFFFFFFFF) >> 1) as u32;
            //let strand0 = y0 & 0x1;
            let hash = mer0.x >> 8;
            let other = shmmr_map.get(&hash).unwrap();
            if other.len() > 10 || other.len() < 2 {
                continue;
            }
            for y1 in other {
                let rid1 = (*y1 >> 32) as u32;
                if rid1 == *rid {
                    continue;
                }
                let pos1 = ((*y1 & 0xFFFFFFFF) >> 1) as u32;
                //let strand1 = *y1 & 0x1;
                matches
                    .entry((*rid, rid1))
                    .or_insert_with(|| vec![])
                    .push((pos0, pos1));
                //log::info!("S {} {} {} {} {} {} {}", rid, pos0, strain0, rid1, pos1, strain1, (pos0 as i32) - (pos1 as i32));
            }
        }
    }

    let rel_path = format!("{}_rel.dat", output_prefix);
    let mut rel_file = File::create(rel_path).unwrap();
    let mut a_to_p = FxHashMap::<u32, u32>::default();

    for ((rid0, rid1), v) in matches {
        let n0 = id2name.get(&rid0).unwrap();
        let n1 = id2name.get(&rid1).unwrap();

        let s0 = shmmr_db.get(n0).unwrap().1.len() as f32;
        let s1 = shmmr_db.get(n1).unwrap().1.len() as f32;
        let c = v.len();
        let r0 = (c as f32) / s0;
        let r1 = (c as f32) / s1;
        writeln!(
            rel_file,
            "S {} {} {} {} {} {} {}",
            rid0, rid1, c, s0, s1, r0, r1
        )?;
        if r0 > 0.50 && s1 > s0 {
            a_to_p.insert(rid0, rid1);
        }
    }

    let p_ctg_path = format!("{}_p.fa", output_prefix);
    let mut p_ctg_file = File::create(p_ctg_path).unwrap();
    let a_ctg_path = format!("{}_a.fa", output_prefix);
    let mut a_ctg_file = File::create(a_ctg_path).unwrap();
    let mut ctg_ids = id2name.keys().map(|x| *x).collect::<Vec<u32>>();
    ctg_ids.sort();
    let mut p_ctg_lengths = Vec::<(String, usize)>::new();
    let mut a_ctg_lengths = Vec::<(String, usize)>::new();
    for ctg_id in ctg_ids {
        if a_to_p.contains_key(&ctg_id) {
            let n0 = id2name.get(&ctg_id).unwrap();
            let n1 = id2name.get(a_to_p.get(&ctg_id).unwrap()).unwrap();
            writeln!(rel_file, "A {} {}", n0, n1)?;
            let seq = seq_db.get(n0).unwrap();
            writeln!(a_ctg_file, ">{}", n0)?;
            writeln!(a_ctg_file, "{}", String::from_utf8_lossy(seq))?;
            a_ctg_lengths.push((n0.clone(), seq.len()))
        } else {
            let n0 = id2name.get(&ctg_id).unwrap();
            writeln!(rel_file, "P {} {}", n0, n0)?;
            let seq = seq_db.get(n0).unwrap();
            writeln!(p_ctg_file, ">{}", n0)?;
            writeln!(p_ctg_file, "{}", String::from_utf8_lossy(seq))?;
            p_ctg_lengths.push((n0.clone(), seq.len()))
        }
    }
    log::info!("primary ctg stats");
    log_asm_summary(p_ctg_lengths);
    log::info!("associated ctg stats");
    log_asm_summary(a_ctg_lengths);
    Ok(())
}
