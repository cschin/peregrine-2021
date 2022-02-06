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
// take the layout file and generate contig sequence from the layout file
// it handles the read stitches and consensuse
//

use super::shmmrutils::match_reads;
use super::shmmrutils::sequence_to_shmmrs;
use super::shmmrutils::{get_2bit_fragment, ReadLocation};
use memmap::MmapOptions;
use std::fs::File;
use std::io::{self, BufRead, Write};
use std::path::Path;

use rustc_hash::{FxHashMap, FxHashSet};

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn get_shmr_offset(s0: &Vec<u8>, s1: &Vec<u8>) -> (u32, u32) {
    
    // take two sequences and use the shmmrs to compute/etimate the offset between them
    // here we take the first hit, it might be useful to use the most common offset 

    let mmer0 = sequence_to_shmmrs(0, &s0, 24, 24, 1);
    let mmer1 = sequence_to_shmmrs(0, &s1, 24, 24, 1);
    let mut counter0 = FxHashMap::<u64, u32>::default();
    for m in mmer0.iter() {
        let x = m.x >> 8;
        *counter0.entry(x).or_insert(0) += 1;
    }
    let mut hmap0 = FxHashMap::<u64, u32>::default();
    for m in mmer0.iter() {
        let x = m.x >> 8;
        let pos = ((m.y & 0xFFFFFFFF) >> 1) as u32;
        if *counter0.get(&x).unwrap() == 1 {
            hmap0.insert(x, pos);
        }
    }

    let (mut offset0, mut offset1) = (0_u32, 0_u32);
    for m in mmer1 {
        let x = m.x >> 8;
        let pos = ((m.y & 0xFFFFFFFF) >> 1) as u32;
        if hmap0.contains_key(&x) {
            offset0 = *hmap0.get(&x).unwrap();
            offset1 = pos;
            break;
        }
    }
    (offset0, offset1)
}

fn get_ctg_cns_with_tiling_reads(
    tiling_reads: &Vec<(u32, Vec<u8>)>,
    out_frag_bgn: &Vec<(usize, usize)>,
    template_seq: &Vec<u8>,
) -> (Vec<u8>, f32) {
    //
    // generate consensus squence from tiling reads
    //
    let mut cov_aux = Vec::<(usize, i32)>::new();
    let mut s_bgns = FxHashSet::<usize>::default();
    let mut s_ends = FxHashSet::<usize>::default();
    let mut delta_count = FxHashMap::<(usize, u8, u8), u32>::default();
    for i in 0..tiling_reads.len() {
        let (_rid0, seq0) = tiling_reads.get(i).unwrap();
        let (bgn, bgn0) = *out_frag_bgn.get(i).unwrap();
        let len = seq0.len() - bgn0;
        let end = if bgn + len < template_seq.len() {
            bgn + len
        } else {
            template_seq.len()
        };
        let end0 = if bgn0 + len < seq0.len() {
            bgn + len
        } else {
            seq0.len()
        };

        let seq_s = template_seq[bgn..end].to_vec();
        let seq0_s = seq0[bgn0..end0].to_vec();
        let ovlp = match_reads(&seq_s, &seq0_s, true, 0.02, 1200, 256);

        if let Some(ovlp) = ovlp {
            //println!("{} {} {} {}", ovlp.bgn0, ovlp.end0, ovlp.bgn1, ovlp.end1);
            let mut dpts = ovlp.deltas.unwrap();
            dpts.reverse();
            let mut d = 0_u8;
            let mut px = 0_usize;
            cov_aux.push((bgn + ovlp.bgn0 as usize, 1));
            cov_aux.push((bgn + ovlp.end0 as usize, -1));
            s_bgns.insert(bgn + ovlp.bgn0 as usize);
            s_ends.insert(bgn + ovlp.end0 as usize);
            for dpt in dpts {
                let c = if dpt.dk > 0 {
                    b'-'
                } else {
                    seq0_s[dpt.y as usize - 1]
                };
                let cx = dpt.x as usize + bgn - 1;
                cov_aux.push((cx, 0));
                if cx != px {
                    d = 0;
                } else {
                    d += 1;
                }
                //println!("D {} {} {}", cx, d, c as char);
                *delta_count.entry((cx, d, c)).or_insert(0) += 1;
                px = cx;
            }
        }
    }
    cov_aux.sort();
    let mut dpt_cov = FxHashMap::<usize, u32>::default();
    let mut cov = 0_i32;
    for (p, d) in cov_aux {
        if d != 0 {
            cov += d;
        } else {
            dpt_cov.entry(p).or_insert(cov as u32);
        }
    }

    let mut dpt_c = Vec::<(usize, u8, u8)>::new();
    for k in delta_count.keys() {
        let count = delta_count.get(k).unwrap();
        let cov = dpt_cov.get(&k.0).unwrap();
        if *count > (cov >> 1) {
            /*
            println!(
                "D {} {} {} {} {} {} {}",
                k.0,
                k.1,
                template_seq[k.0] as char,
                k.2 as char,
                count,
                cov,
                String::from_utf8_lossy(&template_seq[k.0..k.0 + 5].to_vec())
            );
            */
            dpt_c.push(*k);
        }
    }
    dpt_c.sort();

    let mut out_seq = Vec::<u8>::new();
    let mut pos_map = Vec::<usize>::new();
    let mut c_bgn = 0_usize;
    for k in dpt_c {
        if k.0 > c_bgn {
            if c_bgn != 0 {
                c_bgn += 1;
            }
            let c_end = k.0;
            out_seq.extend(template_seq[c_bgn..c_end].to_vec());
            pos_map.extend((c_bgn..c_end).collect::<Vec<usize>>());
            c_bgn = k.0;
        }
        if k.0 == c_bgn {
            if k.1 == 0 && k.2 != b'-' {
                out_seq.push(template_seq[k.0]);
                pos_map.push(k.0);
            }
            if k.2 != b'-' {
                out_seq.push(k.2);
                pos_map.push(k.0);
            }
        }
    }
    let c_end = template_seq.len();
    out_seq.extend(template_seq[c_bgn..c_end].to_vec());
    pos_map.extend((c_bgn..c_end).collect::<Vec<usize>>());
    let mut cov = 0_u32;
    let mut sidx = 0;
    let mut lq_count = 0_u32;
    let out_seq2 = pos_map
        .iter()
        .map(|&p| {
            if s_bgns.contains(&p) {
                cov += 1;
                s_bgns.remove(&p);
            }
            if s_ends.contains(&p) {
                cov -= 1;
                s_ends.remove(&p);
            }
            let base = out_seq[sidx];
            sidx += 1;

            if cov < 3 {
                lq_count += 1;
                base + 4
            } else {
                base
            }
        })
        .collect::<Vec<u8>>();
    /*
        println!(
        "{:?} {:?} {:?}",
        lq_count,
        out_seq2.len(),
        (lq_count as f32 / out_seq2.len() as f32)
    );
    */
    let r = (lq_count as f32) / (out_seq2.len() as f32);
    (out_seq2, r)
}

fn stitch_fragments(
    tiling_reads: &Vec<(u32, Vec<u8>)>,
    match_bng: &FxHashMap<(u32, u32), (u32, u32)>,
) -> (Vec<u8>, f32) {
    let mut frag_bgn = vec![0_u32; tiling_reads.len()];
    let mut frag_end = vec![0_u32; tiling_reads.len()];
    for i in 0..tiling_reads.len() - 1 {
        let (rid0, seq0) = tiling_reads.get(i).unwrap();
        let (rid1, seq1) = tiling_reads.get(i + 1).unwrap();
        //println!("DBG {} {}", rid0, rid1);
        let (mut bgn0, mut bgn1) = match_bng.get(&(*rid0, *rid1)).unwrap();
        if bgn0 < frag_bgn[i] as u32 {
            let correction = frag_bgn[i] as u32 - bgn0;
            bgn0 += correction;
            bgn1 += correction;
        }
        let end0 = if bgn0 as usize + 100 < seq0.len() {
            bgn0 as usize + 100
        } else {
            seq0.len()
        };
        let end1 = if bgn1 as usize + 100 < seq1.len() {
            bgn1 as usize + 100
        } else {
            seq1.len()
        };
        let seq0_s = seq0[bgn0 as usize..end0].to_vec();
        let seq1_s = seq1[bgn1 as usize..end1].to_vec();
        let (offset0, offset1) = get_shmr_offset(&seq0_s, &seq1_s);
        //println!("F {} {} {} {}", rid0, rid1, offset0, offset1);
        frag_end[i] = bgn0 + offset0;
        frag_bgn[i + 1] = bgn1 + offset1;
    }
    frag_end[tiling_reads.len() - 1] =
        tiling_reads.get(tiling_reads.len() - 1).unwrap().1.len() as u32;

    let mut template_seq = Vec::<u8>::new();
    let mut out_frag_bgn = vec![(0_usize, 0_usize); tiling_reads.len()];
    for i in 0..tiling_reads.len() {
        let (_rid0, seq0) = tiling_reads.get(i).unwrap();
        let b = *frag_bgn.get(i).unwrap() as usize;
        let e = *frag_end.get(i).unwrap() as usize;
        out_frag_bgn[i] = (template_seq.len(), b);
        template_seq.extend(seq0[b..e].to_vec());
    }

    get_ctg_cns_with_tiling_reads(&tiling_reads, &out_frag_bgn, &template_seq)
}

pub fn log_asm_summary(ctg_lengths: Vec<(String, usize)>) -> () {
    let mut lengths = ctg_lengths.iter().map(|x| x.1).collect::<Vec<usize>>();
    lengths.sort();
    lengths.reverse();
    let total_bases: usize = lengths.iter().sum();
    let total_ctgs = lengths.len();
    let mut n50 = 0_usize;
    let mut n90 = 0_usize;
    let mut count_gt_100kb = 0_u32;
    log::info!("Total size: {}", total_bases);
    log::info!("Longest size: {}", lengths.get(0).unwrap_or(&0));
    let mut cumsum = 0_usize;
    for l in lengths {
        cumsum += l;
        if cumsum as f32 > (total_bases as f32 * 0.5) {
            if n50 == 0 {
                n50 = l;
            }
        }
        if cumsum as f32 > (total_bases as f32 * 0.9) {
            if n90 == 0 {
                n90 = l;
            }
        }
        if l > 100000 {
            count_gt_100kb += 1;
        }
    }
    log::info!("N50: {}", n50);
    log::info!("N90: {}", n90);
    log::info!("Number of contigs: {}", total_ctgs);
    log::info!("Number of Contigs > 100kb: {}", count_gt_100kb);
}

pub fn layout2ctg(
    seqdb_file: &String,
    index_file: &String,
    layout_file: &String,
    output_file_prefix: &String,
) -> Result<(), io::Error> {
    let mut read_index = Vec::<ReadLocation>::new();

    if let Ok(lines) = read_lines(index_file) {
        for line in lines {
            if let Ok(rec) = line {
                // let rec_trimmed = rec.trim_end();
                // the record line looks like 000000023 m64062_190803_042216/144/ccs 20359 467415
                let v: Vec<&str> = rec.split_whitespace().collect();
                // let rid: u32 = v[0].parse().unwrap();
                let start: usize = v[3].parse().unwrap();
                let len: usize = v[2].parse().unwrap();
                read_index.push(ReadLocation {
                    start: start,
                    len: len,
                });
                // println!("{} {} {}", rid, start, len);
            }
        }
    }

    let file = File::open(seqdb_file).unwrap();
    let readsdb = unsafe { MmapOptions::new().map(&file).unwrap() };

    let mut p_ctg_file = File::create(&format!("{}_m.fa", output_file_prefix)).unwrap();
    let mut a_ctg_file = File::create(&format!("{}_e0.fa", output_file_prefix)).unwrap();
    let mut tiling_reads = Vec::<(u32, Vec<u8>)>::new();
    let mut match_bgn = FxHashMap::<(u32, u32), (u32, u32)>::default();

    let mut pre_ctg_id = 0_u32;
    let mut pre_ctg_tag = 'P';
    let mut ctg_lengths = Vec::<(String, usize)>::new();
    let basemap = [b'A', b'C', b'G', b'T', b'a', b'c', b'g', b't'];
    if let Ok(lines) = read_lines(layout_file) {
        for line in lines {
            if let Ok(rec) = line {
                //let rec_trimmed = rec.trim_end();
                let v: Vec<&str> = rec.split_whitespace().collect();
                if v[0] == "P" || v[0] == "D" || v[0] == "A" {
                    if tiling_reads.len() != 0 {
                        let ctg_id = format!("ctg{:06}_{}", pre_ctg_id, pre_ctg_tag);
                        let (bseq, r) = stitch_fragments(&tiling_reads, &match_bgn);
                        if r < 0.6 {
                            let ctg = bseq
                                .iter()
                                .map(|&t| basemap[t as usize])
                                .collect::<Vec<u8>>();
                            match pre_ctg_tag {
                                'P' => {
                                    let _res = writeln!(p_ctg_file, ">{}", ctg_id);
                                    let _res =
                                        writeln!(p_ctg_file, "{}", String::from_utf8_lossy(&ctg));

                                    ctg_lengths.push((ctg_id, ctg.len()));
                                }
                                'A' => {
                                    let _res = writeln!(a_ctg_file, ">{}", ctg_id);
                                    let _res =
                                        writeln!(a_ctg_file, "{}", String::from_utf8_lossy(&ctg));
                                }
                                'D' => {
                                    let _res = writeln!(a_ctg_file, ">{}", ctg_id);
                                    let _res =
                                        writeln!(a_ctg_file, "{}", String::from_utf8_lossy(&ctg));
                                }
                                _ => (),
                            }
                        }
                    }
                    tiling_reads.clear();
                    pre_ctg_tag = v[0].chars().nth(0).unwrap();
                    continue;
                }
                if v[0] == "E" {
                    let utg_id: u32 = v[1].parse().unwrap();
                    let rid0: u32 = v[2].parse().unwrap();
                    let strand0: u8 = v[3].parse().unwrap();
                    let rid1: u32 = v[4].parse().unwrap();
                    let strand1: u8 = v[5].parse().unwrap();
                    let bgn0: u32 = v[6].parse().unwrap();
                    let bgn1: u32 = v[7].parse().unwrap();
                    if tiling_reads.len() == 0 {
                        let rloc = read_index[rid0 as usize];
                        let len = rloc.len as u32;
                        let seq_full =
                            get_2bit_fragment(rid0, strand0, 0, len, &readsdb, &read_index);
                        tiling_reads.push((rid0, seq_full));
                    }
                    let rloc = read_index[rid1 as usize];
                    let len = rloc.len as u32;
                    //let seq_frag = get_seq_fragment(rid1, strand1, bgn, end, &mmap, &read_index);
                    let seq_full = get_2bit_fragment(rid1, strand1, 0, len, &readsdb, &read_index);
                    pre_ctg_id = utg_id;
                    tiling_reads.push((rid1, seq_full));
                    match_bgn.insert((rid0, rid1), (bgn0, bgn1));
                }
            }
        }

        if tiling_reads.len() != 0 {
            let (bseq, r) = stitch_fragments(&tiling_reads, &match_bgn);
            if r < 0.6 {
                let ctg = bseq
                    .iter()
                    .map(|&t| basemap[t as usize])
                    .collect::<Vec<u8>>();
                let ctg_id = format!("ctg{:06}_{}", pre_ctg_id, pre_ctg_tag);
                match pre_ctg_tag {
                    'P' => {
                        let _res = writeln!(p_ctg_file, ">{}", ctg_id);
                        let _res = writeln!(p_ctg_file, "{}", String::from_utf8_lossy(&ctg));

                        ctg_lengths.push((ctg_id, ctg.len()));
                    }
                    'A' => {
                        let _res = writeln!(a_ctg_file, ">{}", ctg_id);
                        let _res = writeln!(a_ctg_file, "{}", String::from_utf8_lossy(&ctg));
                    }
                    'D' => {
                        let _res = writeln!(a_ctg_file, ">{}", ctg_id);
                        let _res = writeln!(a_ctg_file, "{}", String::from_utf8_lossy(&ctg));
                    }
                    _ => (),
                }
            }
            tiling_reads.clear();
        }
    }
    log_asm_summary(ctg_lengths);
    Ok(())
}
