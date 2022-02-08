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
// using the SHIMMER index to find and verify overalaps
//
use super::shmmrutils::*;
use super::Parameters;
use super::{getrusage, RUSAGE_THREAD};
use core::mem::MaybeUninit;
use glob::glob;
use memmap::{Mmap, MmapOptions};
use rustc_hash::FxHashMap;
use rustc_hash::FxHashSet;
use std::convert::TryInto;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write};
use std::mem::size_of;
use std::path::Path;
use structview::{u64_le, View};
use threadpool::ThreadPool;

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

#[derive(Clone, Copy)]
struct ReadLocation {
    start: usize,
    len: usize,
}

#[derive(Clone, Copy)]
struct HitCluster {
    bgn: i32,
    end: i32,
    count: u32,
    strand: u8,
}

#[derive(Clone, Copy, View)]
#[repr(C)]
struct MM128 {
    x: u64_le,
    y: u64_le,
}

type DeltaMap = FxHashMap<(u32, u8, i32), Vec<u32>>;

impl fmt::Display for MM128 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let x = self.x.to_int();
        let y = self.y.to_int();
        let kmer = x >> 8;
        let span = x & 0xFF;
        let rid = y >> 32;
        let last_pos = (y & 0xFFFFFFFF) >> 1;
        let strand = y & 0x01;
        write!(f, "{} {} {} {} {}", kmer, span, rid, last_pos, strand)
    }
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

fn format_overlap(o: Overlap) -> String {
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

fn get_marker_cov(mut ovlps: Vec<Overlap>, delta_map: &DeltaMap) -> Vec<Overlap> {
    let mut anchors: Vec<u32> = Vec::with_capacity(256);
    let mut delta_count: FxHashMap<u32, u32> = FxHashMap::default();
    delta_count.reserve(256);
    //let rid0 = ovlps[0].rid0;

    ovlps.sort_by(|a, b| b.idt.partial_cmp(&a.idt).unwrap());
    let mut rid1set: FxHashSet<u32> = FxHashSet::default();
    rid1set.reserve(256);
    let mut new_ovlps: Vec<Overlap> = Vec::with_capacity(256);
    for v in ovlps.iter_mut() {
        if !rid1set.contains(&v.rid1) {
            rid1set.insert(v.rid1);
            new_ovlps.push(*v);
        }
    }

    for v in new_ovlps.iter() {
        if v.bgn0 > 5 {
            anchors.push((v.bgn0 - 5) << 2); //note vv.bgn0 is an unsigned type
        } else {
            anchors.push(0);
        }
        anchors.push((v.end0 + 5) << 2 | 0x3);
        let key = &(v.rid1, v.strand1, v.d_left);
        if let Some(deltas) = delta_map.get(&key) {
            // reduce duplicate delta point in single sequence first
            let mut delta_set = FxHashSet::<u32>::default();
            delta_set.reserve(256);
            for d in deltas.iter() {
                delta_set.insert(*d);
            }
            for d in delta_set.iter() {
                anchors.push(*d);
                let c = delta_count.entry(*d).or_insert(0);
                *c += 1;
            }
        }
    }
    anchors.sort_unstable();
    let mut cov: u32 = 0;
    let mut dcount: u32 = 0;
    let mut marker_set = FxHashSet::<u32>::default();
    marker_set.reserve(256);
    for v in anchors.iter() {
        let anchor_type: u32 = v & 0x3;
        match anchor_type {
            0 => {
                cov += 1;
                dcount = 0;
            }
            1 | 2 => {
                dcount = *delta_count.get(v).unwrap();
            }
            3 => {
                cov -= 1;
                dcount = 0;
            }
            _ => (),
        }
        if dcount <= cov && anchor_type != 3 && anchor_type != 0 {
            if dcount > 4 && cov - dcount > 4 {
                marker_set.insert(*v);
            }
        }
    }

    for v in new_ovlps.iter_mut() {
        let mut marker_count: u32 = 0;
        let mut max_marker_count: u32 = 0;
        let key = &(v.rid1, v.strand1, v.d_left);
        if let Some(deltas) = delta_map.get(&key) {
            let mut delta_set = FxHashSet::<u32>::default();
            delta_set.reserve(256);
            for d in deltas.iter() {
                delta_set.insert(*d);
            }
            for d in delta_set.iter() {
                if marker_set.contains(d) {
                    marker_count += 1;
                }
            }
            for d in marker_set.iter() {
                let p = d >> 2;
                if p > v.bgn0 && p <= v.end0 {
                    max_marker_count += 1;
                }
            }
        }
        (*v).dist_c = marker_count;
        (*v).max_dist_c = max_marker_count;
        (*v).idt_c = 100.0 - 100.0 * (marker_count as f32) / ((v.end0 - v.bgn0) as f32);
    }

    new_ovlps.sort_by(|a, b| a.dist_c.cmp(&b.dist_c));

    let mut rp = FxHashSet::<u32>::default();
    rp.reserve(256);
    for v in new_ovlps.iter_mut() {
        if rp.contains(&v.rid1) {
            continue;
        }
        rp.insert(v.rid1);
        if (v.max_dist_c > 2 && v.dist_c <= 1) || (v.max_dist_c <= 2 && v.dist_c == 0) {
            v.flag |= 0x2;
        }
    }

    new_ovlps
}

fn find_offset(s0: &Vec<u8>, s1: &Vec<u8>, s0_range: (u32, u32)) -> Option<(u32, u32)> {
    // find offset using 8-kmer match
    let mut s0map = FxHashMap::<u32, u32>::default(); //hash to position
    s0map.reserve(256);
    let mut h = 0_u32;
    let bgn = if (s0_range.0 as i32) - 100 > 0 {
        s0_range.0 - 100
    } else {
        0
    };
    let end = if s0_range.1 + 100 < s0.len() as u32 {
        s0_range.1 + 100 - bgn
    } else {
        s0.len() as u32 - bgn
    };
    //println!("s0_range {} {}", s0_range.0, s0_range.1);
    //println!("bgn/end {} {}", bgn, end);
    for i in 0..8 {
        h = (h << 4) | (s0[(bgn + i) as usize] & 0b0011) as u32;
    }
    s0map.insert(h, bgn + 1);
    for i in 8..end {
        h = (h << 4) | (s0[(bgn + i) as usize] & 0b0011) as u32;
        s0map.entry(h).or_insert(bgn + i - 8 + 1);
    }

    h = 0_u32;
    for i in 0..8 {
        h = (h << 4) | (s1[i as usize] & 0b0011) as u32;
    }
    let end = 48;
    let mut offset = 0_u32;
    let mut offset2 = 0_u32;
    let mut found = false;
    for i in 8..end {
        if s0map.contains_key(&h) {
            offset = *s0map.get(&h).unwrap();
            offset2 = i - 8;
            found = true;
            break;
        }
        h = (h << 4) | (s1[i as usize] & 0b0011) as u32;
    }
    if found {
        Some((offset, offset2))
    } else {
        None
    }
}

fn process_hits(
    rid0: u32,
    rid1: u32,
    c: &HitCluster,
    read_index: &Vec<ReadLocation>,
    readsdb: &Mmap,
    tol: f64,
) -> Option<(Overlap, Option<Vec<DeltaPoint>>)> {
    let strand: u8 = c.strand;
    let rloc0 = read_index[rid0 as usize];
    let rloc1 = read_index[rid1 as usize];
    let mut bgn = c.bgn;
    let mut end = c.end;
    let s0 = rloc0.start;
    let len0 = rloc0.len;
    let e0 = s0 + len0;
    let s1 = rloc1.start;
    let len1 = rloc1.len;
    let e1 = s1 + len1;
    // println!("DBG: {} {} {} {}", rid0, rid1, len0, len1) ;
    let mut seq0 = Vec::<u8>::with_capacity(len0);
    let mut seq1 = Vec::<u8>::with_capacity(len1);

    for c in &readsdb[s0..e0] {
        let base_flag_2bit = c & 0x0F;
        seq0.push(base_flag_2bit);
    }

    if strand == 0 {
        for c in &readsdb[s1..e1] {
            let base_flag_2bit = c & 0x0F;
            seq1.push(base_flag_2bit);
        }
    } else {
        for c in &readsdb[s1..e1] {
            let base_flag_2bit = (c >> 4) & 0x0F;
            seq1.push(base_flag_2bit);
        }
        bgn = bgn - (rloc1.len as i32);
        end = end - (rloc1.len as i32);
    }

    let offset0: u32;
    let offset1: u32;

    if bgn >= 0 {
        if let Some(tmp) = find_offset(&seq0, &seq1, (bgn as u32, end as u32)) {
            //println!("offset {} {} {} {} {} {}", rid0, rid1, bgn, end, tmp.0, tmp.1);
            offset0 = tmp.0;
            offset1 = tmp.1;
        } else {
            return None;
        }
    } else {
        assert!(-end <= -bgn);
        //assert!(-end > 0, "{}", end);
        //assert!(-bgn > 0, "{}", bgn);
        if end > 0 {
            end = 0
        };

        if let Some(tmp) = find_offset(&seq1, &seq0, (-end as u32, -bgn as u32)) {
            //println!("offset {} {} {} {} {} {}", rid0, rid1, bgn, end, tmp.0, tmp.1);
            offset1 = tmp.0;
            offset0 = tmp.1;
        } else {
            return None;
        }
        //let tmp = find_offset(&seq1, &seq0, (-bgn as u32, -end as u32));
        //println!("offset {} {} {} {}", bgn, end, offset0, tmp);
    }

    //assert!((offset1 as usize) < seq1.len(), "DDD {} {} {} {} {} {}", rid0, rid1, len0, len1, bgn, end);
    let seq0 = get_hpc_seq(&seq0[(offset0 as usize)..].to_vec());
    let seq1 = get_hpc_seq(&seq1[(offset1 as usize)..].to_vec());

    if let Some(ovlpmatch) = match_reads(&seq0.s, &seq1.s, true, tol, 1200, 16) {
        let idt_est = 100.0_f32 - 100.0 * (ovlpmatch.dist as f32) / (ovlpmatch.m_size as f32);
        let left = bgn;
        let right = bgn + (len1 as i32) - (len0 as i32);
        let overlap = Overlap {
            rid0: rid0,
            rid1: rid1,
            strand1: strand,
            len0: len0 as u32,
            len1: len1 as u32,
            d_left: left,
            d_right: right,
            bgn0: seq0.p[ovlpmatch.bgn0 as usize] + offset0,
            end0: seq0.p[ovlpmatch.end0 as usize] + offset0,
            bgn1: seq1.p[ovlpmatch.bgn1 as usize] + offset1,
            end1: seq1.p[ovlpmatch.end1 as usize] + offset1,
            dist: ovlpmatch.dist,
            idt: idt_est,
            dist_c: 0_u32,
            max_dist_c: 0_u32,
            idt_c: 0_f32,
            flag: 0x00,
        };
        let mut new_deltas = Vec::<DeltaPoint>::with_capacity(256);
        for d in ovlpmatch.deltas.unwrap().iter() {
            let x = seq0.p[d.x as usize] + offset0;
            let y = seq1.p[d.y as usize] + offset1;
            let dk = d.dk;
            if x > overlap.bgn0 && x < overlap.end0 {
                new_deltas.push(DeltaPoint { x: x, y: y, dk: dk });
            }
        }
        Some((overlap, Some(new_deltas)))
    } else {
        None
    }
}

fn extract_candidates(
    rid: &u32,
    shmmrs: &Vec<MM128>,
    smp_index: &FxHashMap<u128, Vec<u64>>,
) -> (FxHashMap<u32, Vec<(i32, bool)>>, u32, u32) {
    let mmer0 = shmmrs.get(0).unwrap();
    let x0 = mmer0.x.to_int();
    let mut y0 = mmer0.y.to_int();
    //let mut rid0 = (y0 >> 32) as u32;
    let mut hash0 = (x0 >> 8) as u128;
    let mut candidates = FxHashMap::<u32, Vec<(i32, bool)>>::default();
    candidates.reserve(256);
    let smp_count_limit = 12000;
    let mut overlimit_count = 0;
    let mut total_count = 0;
    for i in 1..shmmrs.len() {
        let mmer1 = shmmrs.get(i).unwrap();
        let x1 = mmer1.x.to_int();
        let y1 = mmer1.y.to_int();
        //let rid1 = (y1 >> 32) as u32;
        let hash1 = (x1 >> 8) as u128;

        if hash0 != hash1 {
            //forward
            let h128 = (hash0 << 64) | hash1;
            if let Some(counter) = smp_index.get(&h128) {
                total_count += 1;
                let p0 = ((y0 & 0xFFFFFFFF) >> 1) as i32;
                //let strand0 = y0 & 0x01;
                //let span = (x1 & 0xFF) as i32;
                if counter.len() < smp_count_limit {
                    for m_rid in counter.iter() {
                        let rid1 = (m_rid >> 32) as u32;
                        if rid1 == *rid {
                            continue;
                        }
                        let p = ((m_rid & 0xFFFFFFFF) >> 1) as i32;
                        let reversed = false;
                        // TODO (July 25, 2021): check the mmer strand and pair strand is consistent
                        // We don't run into problem as we check the alignment later, maybe it is beneficial to fitering those in-consistent one firts
                        // if  (m_rid & 0x01) != (y0 & 0x01) { continue };
                        candidates
                            .entry(rid1)
                            .or_insert_with(|| Vec::<(i32, bool)>::with_capacity(1024))
                            .push((p0 - p, reversed));
                        //println!("1: {} {} {} {} {} {} {}", rid, p0, rid1, p0 - p, strand0, (m_rid & 0x01), reversed);
                    }
                } else {
                    overlimit_count += 1;
                }
            }
            //reverse
            let h128 = (hash1 << 64) | hash0;
            if let Some(counter) = smp_index.get(&h128) {
                total_count += 1;
                let p0 = ((y1 & 0xFFFFFFFF) >> 1) as i32;
                //let strand0 = y1 & 0x01;
                let span = (x1 & 0xFF) as i32;
                if counter.len() < smp_count_limit {
                    for m_rid in counter.iter() {
                        let rid1 = (m_rid >> 32) as u32;
                        if rid1 == *rid {
                            continue;
                        }
                        let p = ((m_rid & 0xFFFFFFFF) >> 1) as i32;
                        let reversed = true;
                        // TODO (July 25, 2021): check the mmer strand and pair strand is consistent
                        // We don't run into problem as we check the alignment later, maybe it is beneficial to fitering those in-consistent one firts
                        // if  (m_rid & 0x01) == (y0 & 0x01) { continue };
                        candidates
                            .entry(rid1)
                            .or_insert_with(|| Vec::<(i32, bool)>::with_capacity(1024))
                            .push((p0 + p - span + 1, reversed));
                        //println!("2: {} {} {} {} {} {} {}", rid, p0, rid1, p + p0 - span + 1, strand0, (m_rid & 0x01), reversed);
                    }
                } else {
                    overlimit_count += 1;
                }
            }
        }
        //x0 = x1;
        y0 = y1;
        //rid0 = rid1;
        hash0 = hash1;
    }
    (candidates, total_count, overlimit_count)
}

fn extract_hit_cluster(hits: &mut Vec<(i32, bool)>) -> Vec<HitCluster> {
    let mut clusters = Vec::<HitCluster>::with_capacity(256);

    let mut cluster_end = -1000000_i32;
    let mut cluster_start = -1000000_i32;
    let mut item_count = 0_u32;
    let item_count_limit = 2;
    hits.sort_by(|a, b| a.0.cmp(&b.0));
    for (pt, rev) in hits.iter() {
        if *rev == true {
            continue;
        }
        //println!("Y0 {} {} {} {} {} {}", rid, rid1, pt, cluster_end, *pt-cluster_end, rev);
        item_count += 1;
        if pt - cluster_end > 24 {
            if item_count > item_count_limit {
                clusters.push(HitCluster {
                    bgn: cluster_start,
                    end: cluster_end,
                    count: item_count,
                    strand: 0,
                });
                //println!("C {} {} {} {} {} {}", rid, rid1, cluster_start, cluster_end, item_count, 0);
            }
            cluster_start = *pt;
            item_count = 0;
        }
        cluster_end = *pt;
    }
    if item_count > item_count_limit {
        clusters.push(HitCluster {
            bgn: cluster_start,
            end: cluster_end,
            count: item_count,
            strand: 0,
        });
        //println!("C {} {} {} {} {} {}", rid, rid1, cluster_start, cluster_end, item_count, 0);
    }

    let mut cluster_end = -1000000_i32;
    let mut cluster_start = -1000000_i32;
    let mut item_count = 0_u32;
    for (pt, rev) in hits.iter() {
        if *rev == false {
            continue;
        }
        //println!("Y0 {} {} {} {} {} {}", rid, rid1, pt, cluster_end, *pt-cluster_end, rev);
        item_count += 1;
        if pt - cluster_end > 24 {
            if item_count > item_count_limit {
                clusters.push(HitCluster {
                    bgn: cluster_start,
                    end: cluster_end,
                    count: item_count,
                    strand: 1,
                });
                //println!("C {} {} {} {} {} {}", rid, rid1, cluster_start, cluster_end, item_count, 1);
            }
            cluster_start = *pt;
            item_count = 0;
        }
        cluster_end = *pt;
    }
    if item_count > item_count_limit {
        clusters.push(HitCluster {
            bgn: cluster_start,
            end: cluster_end,
            count: item_count,
            strand: 1,
        });
        //println!("C {} {} {} {} {} {}", rid, rid1, cluster_start, cluster_end, item_count, 1);
    }
    //println!("X1 {} {} {}", rid, rid1, hits.len());
    clusters
}

fn output_ovlp_candidate_for_chunk(
    smp_index: &FxHashMap<u128, Vec<u64>>,
    read_shmmr_index: &FxHashMap<u32, Vec<MM128>>,
    readsdb: &Mmap,
    read_index: &Vec<ReadLocation>,
    output_prefix: &String,
    total_chunk: u32,
    mychunk: u32,
    tol: f64,
) {
    let min_hit_limit = 2;

    let mut file =
        BufWriter::new(File::create(format!("{}.{:02}", output_prefix, mychunk)).unwrap());

    let mut rdata = unsafe { MaybeUninit::uninit().assume_init() };
    let _res = unsafe { getrusage(RUSAGE_THREAD, &mut rdata) };
    let mut current_utime = rdata.ru_utime.tv_sec;
    let mut current_stime = rdata.ru_stime.tv_sec;
    let mut count = 0_u64;

    for (rid0, shmmrs) in read_shmmr_index.iter() {
        if rid0 % total_chunk != mychunk % total_chunk {
            continue;
        }
        count += 1;
        if count % 10000 == 0 {
            let _res = unsafe { getrusage(RUSAGE_THREAD, &mut rdata) };
            let current_utime_now = rdata.ru_utime.tv_sec;
            let current_stime_now = rdata.ru_stime.tv_sec;
            let dutime = current_utime_now - current_utime;
            let dstime = current_stime_now - current_stime;
            log::info!(
                "ovlp chunk {}, utime: {} s / 10000 reads, total: {} reads",
                mychunk,
                dutime,
                count
            );
            log::info!(
                "ovlp chunk {}, stime: {} s / 10000 reads, total: {} reads",
                mychunk,
                dstime,
                count
            );
            if dstime as f64 / (dutime + dstime) as f64 > 0.3 {
                log::warn!(
                    "excessive system cpu usage, consider to reduce conccurance or incerease RAM"
                );
            }
            current_utime = current_utime_now;
            current_stime = current_stime_now;
        }
        let (mut candidates, total_count, overlimit_count) =
            extract_candidates(rid0, shmmrs, smp_index);
        let _res = writeln!(file, "R0 {} {} {}", rid0, total_count, overlimit_count);
        let total_number_of_cluster: usize;
        let mut tries: u32;
        if candidates.len() > 0 {
            let mut ovlps = Vec::<Overlap>::with_capacity(128);
            let mut delta_map = DeltaMap::default();
            delta_map.reserve(256);

            let mut cluster_v = Vec::<(u32, HitCluster)>::with_capacity(512);
            for (rid1, hits) in candidates.iter_mut() {
                if hits.len() > min_hit_limit {
                    let clusters = extract_hit_cluster(hits);
                    for c in clusters.iter() {
                        cluster_v.push((*rid1, *c));
                    }
                }
            }
            cluster_v.sort_by(|a, b| b.1.count.cmp(&a.1.count));
            total_number_of_cluster = cluster_v.len();

            tries = 0;
            for (rid1, c) in cluster_v.iter() {
                //writeln!(file, "C {} {} {} {} {} {}", rid0, rid1, c.bgn, c.end, c.count, c.strand);
                if let Some((overlap, Some(dpts))) =
                    process_hits(*rid0, *rid1, c, read_index, readsdb, tol)
                {
                    ovlps.push(overlap);
                    let mut new_dpts = Vec::<u32>::with_capacity(128);
                    for dpt in dpts {
                        let d = dpt.x << 2 | (if dpt.dk > 0 { 1 } else { 2 });
                        new_dpts.push(d);
                    }
                    delta_map.insert((*rid1, overlap.strand1, overlap.d_left), new_dpts);
                    //writeln!(file, "X {}", format_overlap(overlap));
                }
                tries += 1;
                if tries > 320 {
                    break;
                }
            }
            let ovlps_len = ovlps.len();

            if ovlps.len() > 0 {
                let mut updated_ovlps = get_marker_cov(ovlps, &delta_map);
                updated_ovlps.sort_by(|a, b| a.d_left.cmp(&b.d_left));
                for overlap in updated_ovlps {
                    let _res = writeln!(file, "O {}", format_overlap(overlap));
                }
            }
            let _res = writeln!(
                file,
                "R1 {} {} {} {}",
                rid0, total_number_of_cluster, tries, ovlps_len
            );
        }
        let _res = writeln!(file, "-");
    }
}

struct ShmrIdxSet {
    filepath: String,
    n_mmers: usize,
    idx_mmap: Mmap,
    cur_mmer_idx: usize,
}

impl ShmrIdxSet {
    fn new(path: &String) -> Result<ShmrIdxSet, io::Error> {
        let file = File::open(path)?;
        let idx_mmap = unsafe { MmapOptions::new().map(&file)? };

        let us = size_of::<usize>();
        let mm128s = size_of::<MM128>();
        let n_mmers = usize::from_le_bytes(idx_mmap[0..us].try_into().unwrap());
        assert!(us + n_mmers * mm128s == idx_mmap.len());

        let cur_mmer_idx = 0;
        Ok(ShmrIdxSet {
            filepath: path.clone(),
            n_mmers: n_mmers,
            idx_mmap: idx_mmap,
            cur_mmer_idx: cur_mmer_idx,
        })
    }
}

impl Iterator for ShmrIdxSet {
    type Item = (u32, Vec<MM128>);
    fn next(&mut self) -> Option<(u32, Vec<MM128>)> {
        let us = size_of::<usize>();
        let mm128s = size_of::<MM128>();

        let mut c_mmer_idx = self.cur_mmer_idx;
        if c_mmer_idx >= self.n_mmers {
            return None;
        }
        let mut out = Vec::<MM128>::with_capacity(256);
        let s = us + mm128s * c_mmer_idx;
        let e = s + mm128s;
        let mmer = MM128::view(&self.idx_mmap[s..e]).unwrap();
        let y0 = mmer.y.to_int();
        let rid = (y0 >> 32) as u32;
        out.push(*mmer);
        c_mmer_idx += 1;

        loop {
            if c_mmer_idx >= self.n_mmers {
                break;
            }
            let s = us + mm128s * c_mmer_idx;
            let e = s + mm128s;
            let mmer = MM128::view(&self.idx_mmap[s..e]).unwrap();
            let y0 = mmer.y.to_int();
            let rid1 = (y0 >> 32) as u32;

            if rid1 != rid {
                break;
            }
            out.push(*mmer);
            c_mmer_idx += 1;
        }
        self.cur_mmer_idx = c_mmer_idx;
        Some((rid, out))
    }
}

pub fn ovlp(
    seqdb_file: &String,
    index_file: &String,
    shmmr_index_file_prefix: &String,
    out_prefix: &String,
    parameters: &Parameters,
) -> Result<(), io::Error> {
    let filename = format!("{}-??-of-??.dat", shmmr_index_file_prefix);
    //let filename = "/wd/resolve_overlap_test/asm_ext_k56/1-index/shmr-L2-??-of-20.dat";
    let mut smp_index = FxHashMap::<u128, Vec<u64>>::default();
    smp_index.reserve(65536);
    let mut read_shmmr_index = FxHashMap::<u32, Vec<MM128>>::default();
    read_shmmr_index.reserve(65536);

    for entry in glob(&filename).expect("Failed to read glob pattern") {
        match entry {
            Ok(path) => {
                log::info!("open index file: {}", path.to_string_lossy());
                let shmmr_idx_set = ShmrIdxSet::new(&path.to_string_lossy().to_string())?;
                let _res = log::info!(
                    "reading index file: {} size: {}",
                    path.as_path().display(),
                    shmmr_idx_set.n_mmers
                );

                for (rid, mmers) in shmmr_idx_set {
                    if mmers.len() < 2 {
                        continue;
                    }
                    let mmer0 = mmers[0];
                    let x0 = mmer0.x.to_int();
                    let mut y0 = mmer0.y.to_int();
                    let mut hash0 = (x0 >> 8) as u128;
                    for i in 1..mmers.len() {
                        let mmer1 = mmers[i];
                        let x1 = mmer1.x.to_int();
                        let y1 = mmer1.y.to_int();
                        let hash1 = (x1 >> 8) as u128;
                        let h128 = (hash0 << 64) | hash1;
                        smp_index.entry(h128).or_insert_with(|| vec![]).push(y0);
                        hash0 = hash1;
                        y0 = y1;
                    }
                    read_shmmr_index.insert(rid, mmers);
                }
            }
            Err(e) => log::error!("{:?}", e),
        }
    }
    log::info!("read_shmmr_index size: {}", read_shmmr_index.len());
    log::info!("smp_index: {}", smp_index.len());

    let mut singletons = Vec::<u128>::with_capacity(65536);
    for (h, v) in smp_index.iter() {
        if v.len() < 2 {
            singletons.push(*h);
        }
    }
    log::info!("remove {} singletons from smp_index", singletons.len());
    for h in singletons {
        smp_index.remove(&h);
    }
    smp_index.shrink_to_fit();
    log::info!("smp_index size updated: {}", smp_index.len());

    let mut read_index = Vec::<ReadLocation>::with_capacity(65535);

    log::info!("read seq idx: {}", seqdb_file);
    let lines = read_lines(index_file)?;
    for line in lines {
        if let Ok(rec) = line {
            //let rec_trimmed = rec.trim_end();
            // the record line looks like 000000023 m64062_190803_042216/144/ccs 20359 467415
            let v: Vec<&str> = rec.split_whitespace().collect();
            //let rid: u32 = v[0].parse().unwrap();
            let start: usize = v[3].parse().unwrap();
            let len: usize = v[2].parse().unwrap();
            read_index.push(ReadLocation {
                start: start,
                len: len,
            });
            //println!("{} {} {}", rid, start, len);
        }
    }

    let read_index = read_index;
    log::info!("read seq db: {}", seqdb_file);
    let file = File::open(seqdb_file)?;
    let readsdb = unsafe { MmapOptions::new().map(&file)? };
    let smp_index = std::sync::Arc::new(smp_index);
    let read_shmmr_index = std::sync::Arc::new(read_shmmr_index);
    let readsdb = std::sync::Arc::new(readsdb);
    let read_index = std::sync::Arc::new(read_index);
    let pool = ThreadPool::new(parameters.nthreads as usize);

    for i in 0..parameters.nchunks {
        let smp_index = smp_index.clone();
        let read_shmmr_index = read_shmmr_index.clone();
        let readsdb = readsdb.clone();
        let read_index = read_index.clone();
        let out_prefix = out_prefix.clone();
        let nchunks = parameters.nchunks;
        let tol = parameters.tol;

        pool.execute(move || {
            output_ovlp_candidate_for_chunk(
                &smp_index,
                &read_shmmr_index,
                &readsdb,
                &read_index,
                &out_prefix,
                nchunks,
                i + 1,
                tol,
            );
        });
    }
    pool.join();
    Ok(())
}
