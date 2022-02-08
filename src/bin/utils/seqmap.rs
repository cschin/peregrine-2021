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
// for handling the read mapping
//

use lazy_static::lazy_static;

use super::build_sdb::{encode_biseq, FastxReader};
use super::shmmrutils::MM128;
use super::shmmrutils::{get_hpc_seq, sequence_to_shmmrs};
use core::ops::Range;
use intervaltree::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
pub type MapIntervalRecord = [u32; 8];
pub type MapIntervals = FxHashMap<u32, IntervalTree<u32, MapIntervalRecord>>;
pub type Shmmrs = Vec<Vec<MM128>>;
use flate2::bufread::GzDecoder;
use petgraph::graphmap::DiGraphMap;
use petgraph::unionfind::UnionFind;
use rayon::prelude::*;
use regex::bytes::Regex;
pub struct SeqDB {
    filepath: String,
    seqs: Vec<Vec<u8>>,
    masks: Vec<Vec<(usize, usize)>>,
    shmmrs: Shmmrs,
    seqlen: FxHashMap<String, usize>,
    id2seqname: FxHashMap<u32, String>,
    seqname2id: FxHashMap<String, u32>,
    sequences_loaded: bool,
    shmmrs_built: bool,
}

impl SeqDB {
    pub fn new(filepath: String) -> Self {
        let seqs = Vec::<Vec<u8>>::default();
        let masks = Vec::<Vec<(usize, usize)>>::default();
        let shimmers_db = Vec::<Vec<MM128>>::new();
        let seqlen = FxHashMap::<String, usize>::default();
        let seqname2id = FxHashMap::<String, u32>::default();
        let id2seqname = FxHashMap::<u32, String>::default();

        SeqDB {
            filepath: filepath,
            seqs: seqs,
            masks: masks,
            shmmrs: shimmers_db,
            seqlen: seqlen,
            id2seqname: id2seqname,
            seqname2id: seqname2id,
            sequences_loaded: false,
            shmmrs_built: false,
        }
    }

    pub fn load_sequences(&mut self) -> Result<(), std::io::Error> {
        lazy_static! {
            static ref QR: Regex = Regex::new("[ACGT]+").unwrap();
        }

        let file = File::open(&self.filepath)?;

        let mut reader = BufReader::new(file);
        let mut is_gzfile = false;
        {
            let r = reader.by_ref();
            let mut buf = Vec::<u8>::new();
            let _ = r.take(2).read_to_end(&mut buf);
            if buf == [0x1F_u8, 0x8B_u8] {
                log::info!(
                    "input file: {} detected as gz-compressed file",
                    self.filepath
                );
                is_gzfile = true;
            }
        }

        reader.seek(SeekFrom::Start(0))?;
        if is_gzfile {
            let fastx_buf = BufReader::new(GzDecoder::new(&mut reader));
            let mut fastx_reader = FastxReader::new(fastx_buf, &self.filepath)?;
            let mut sid = 0;
            while let Some(rec) = fastx_reader.next_rec() {
                let rec = rec.unwrap();
                let seqname = String::from_utf8_lossy(&rec.id).into_owned();
                self.id2seqname.insert(sid, seqname.clone());
                self.seqname2id.insert(seqname.clone(), sid);
                self.seqlen.insert(seqname.clone(), rec.seq.len());
                self.seqs.push(encode_biseq(&rec.seq));

                let mask = QR
                    .find_iter(&rec.seq)
                    .map(|m| (m.start(), m.end()))
                    .collect::<Vec<(usize, usize)>>();
                self.masks.push(mask);
                //println!("{:?}", self.mask);
                sid += 1;
            }
        } else {
            let mut fastx_reader = FastxReader::new(reader, &self.filepath).unwrap();
            let mut sid = 0;
            // unfortunatly we, need to repeat the code here as the type fastx_reader is different from the above
            while let Some(rec) = fastx_reader.next_rec() {
                let rec = rec.unwrap();
                let seqname = String::from_utf8_lossy(&rec.id).into_owned();
                self.id2seqname.insert(sid, seqname.clone());
                self.seqname2id.insert(seqname.clone(), sid);
                self.seqlen.insert(seqname.clone(), rec.seq.len());
                self.seqs.push(encode_biseq(&rec.seq));

                let mask = QR
                    .find_iter(&rec.seq)
                    .map(|m| (m.start(), m.end()))
                    .collect::<Vec<(usize, usize)>>();
                self.masks.push(mask);
                //printlnmask!("{:?}", self.mask);
                sid += 1;
            }
        }
        self.sequences_loaded = true;
        Ok(())
    }

    pub fn build_shmmrs(&mut self, w: u32, k: u32, r: u32) -> () {
        let e_seqs = self
            .seqs
            .iter()
            .enumerate()
            .collect::<Vec<(usize, &Vec<u8>)>>();

        let mut out = e_seqs
            .par_chunks(1000)
            .into_par_iter()
            .flat_map(|xv| {
                let mut out = Vec::<(usize, Vec<MM128>)>::new();
                for x in xv.iter() {
                    out.push((x.0, sequence_to_shmmrs(x.0 as u32, x.1, w, k, r)));
                }
                out
            })
            .collect::<Vec<(usize, Vec<MM128>)>>();

        out.par_sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        self.shmmrs.clear();
        for (_sid, shmmrs) in out {
            self.shmmrs.push(shmmrs);
        }

        self.shmmrs_built = true;
        ()
    }

    pub fn get_subseq_by_name(&self, seqname: String, s: usize, e: usize) -> Option<String> {
        if let Some(sid) = self.seqname2id.get(&seqname) {
            if e > *self.seqlen.get(&seqname).unwrap() {
                return None;
            } else {
                let str = String::from_utf8_lossy(&self.seqs.get(*sid as usize).unwrap()[s..e])
                    .to_string();
                Some(str)
            }
        } else {
            None
        }
    }

    pub fn get_subseq_by_id(&self, sid: u32, s: usize, e: usize) -> Option<String> {
        if sid < self.seqs.len() as u32 {
            let seqname = self.id2seqname.get(&sid).unwrap();
            if e > *self.seqlen.get(seqname).unwrap() {
                None
            } else {
                let str = String::from_utf8_lossy(&self.seqs[sid as usize][s..e]).to_string();
                Some(str)
            }
        } else {
            None
        }
    }
}

pub fn load_shmmr_map(mi: &mut MapIntervals, filename: &String) -> Result<(), std::io::Error> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut itvl_vecs = FxHashMap::<u32, Vec<(Range<u32>, MapIntervalRecord)>>::default();

    for r in reader.lines() {
        let r = r?;
        let r: Vec<&str> = r.split_whitespace().collect();
        if r[0] != "M" {
            continue;
        }
        let mut v: MapIntervalRecord = [0; 8];
        for i in 1..v.len() + 1 {
            v[i - 1] = r[i].parse::<u32>().unwrap();
        }
        itvl_vecs
            .entry(v[0])
            .or_insert_with(|| vec![])
            .push((v[1]..v[2], v.clone()));
    }
    for (sid, itvl_vec) in itvl_vecs {
        let itvl: IntervalTree<u32, MapIntervalRecord> = itvl_vec.iter().cloned().collect();
        mi.insert(sid, itvl);
    }
    Ok(())
}

pub fn write_shmmr_map(mi: &MapIntervals, filename: &String) -> Result<(), std::io::Error> {
    let mut file = BufWriter::new(File::create(filename)?);

    for (_rid, itvls) in mi.iter() {
        for itvl in itvls.iter_sorted() {
            let v = itvl.clone().value;
            writeln!(
                file,
                "M {} {} {} {} {} {} {} {}",
                v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]
            )?;
        }
    }
    Ok(())
}

pub fn generate_shmmr_map(shmmrs0: &Shmmrs, shmmrs1: &Shmmrs, max_hits: usize) -> MapIntervals {
    let out_vec = shmmrs1
        .par_chunks(250)
        .into_par_iter()
        .flat_map_iter(|mmers_v| {
            mmers_v
                .iter()
                .filter(|mmers| mmers.len() >= 2)
                .flat_map(|mmers| {
                    let mlen = mmers.len();
                    mmers[0..mlen - 1]
                        .iter()
                        .zip(mmers[1..mlen].iter())
                        .filter(|m| (m.0.x >> 8) != (m.1.x >> 8))
                        .map(|m| {
                            let (m0, m1) = m;
                            let hash0 = (m0.x >> 8) as u128;
                            let hash1 = (m1.x >> 8) as u128;
                            let h128 = (hash0 << 64) | hash1;
                            (h128, m0.y, m1.y)
                        })
                })
        })
        .collect::<Vec<(u128, u64, u64)>>();

    const INDEX_CHUNKS: usize = 16;

    let mut smp_index = Vec::<(u32, FxHashMap<u128, Vec<(u64, u64)>>)>::new();
    for i in 0..INDEX_CHUNKS {
        smp_index.push((i as u32, FxHashMap::<u128, Vec<(u64, u64)>>::default()));
    }
    smp_index.par_iter_mut().for_each(|x| {
        let (idx, idxhash) = x;
        out_vec
            .iter()
            .filter(|&h| ((h.0 & 0xFF) % INDEX_CHUNKS as u128) as u8 == *idx as u8)
            .for_each(|v| {
                idxhash
                    .entry(v.0)
                    .or_insert_with(|| vec![])
                    .push((v.1, v.2))
            })
    });

    let mut all_itvl = MapIntervals::default();

    for (rid0, mmers) in shmmrs0.iter().enumerate() {
        if mmers.len() < 2 {
            continue;
        }

        let mut itvl_vec = Vec::<(Range<u32>, MapIntervalRecord)>::new();

        let mlen = mmers.len();
        if mlen < 2 {
            continue;
        }
        let out_v = mmers[0..mlen - 1]
            .par_iter()
            .zip(mmers[1..mlen].par_iter())
            .filter(|m| (m.0.x >> 8) != (m.1.x >> 8))
            .map(|m| -> [(u128, u64, u64, u32); 2] {
                let (m0, m1) = m;
                let hash0 = (m0.x >> 8) as u128;
                let hash1 = (m1.x >> 8) as u128;
                let h128 = (hash0 << 64) | hash1;
                let rh128 = (hash1 << 64) | hash0;
                [(h128, m0.y, m1.y, 0), (rh128, m0.y, m1.y, 1)]
            })
            .collect::<Vec<[(u128, u64, u64, u32); 2]>>();

        //.collect_into_vec(&mut out_v);

        let mut range_matches = Vec::<Vec<(Range<u32>, MapIntervalRecord)>>::new();
        out_v
            .par_iter()
            .map(|pv| {
                let mut rmv = Vec::<(Range<u32>, MapIntervalRecord)>::new();
                for v in pv {
                    let (h128, y0, y1, orientation) = *v;
                    let shard = ((h128 & 0x0F) % INDEX_CHUNKS as u128) as usize;
                    if !smp_index[shard].1.contains_key(&h128) {
                        continue;
                    }
                    let matches = smp_index[shard].1.get(&h128).unwrap();
                    if (*matches).len() > max_hits {
                        continue;
                    }
                    for (yy0, yy1) in matches.iter() {
                        let rid0 = (y0 >> 32) as u32;
                        let pos0 = ((y0 & 0xFFFFFFFF) >> 1) as u32;
                        let pos1 = ((y1 & 0xFFFFFFFF) >> 1) as u32;

                        let rid1 = (*yy0 >> 32) as u32;
                        let ppos0 = ((*yy0 & 0xFFFFFFFF) >> 1) as u32;
                        let ppos1 = ((*yy1 & 0xFFFFFFFF) >> 1) as u32;
                        let group = 0;

                        match orientation {
                            0 => {
                                if (y0 & 0b01) != (yy0 & 0b01) || (y1 & 0b01) != (yy1 & 0b01) {
                                    continue;
                                }
                            }
                            1 => {
                                if (y0 & 0b01) == (yy1 & 0b01) || (y1 & 0b01) == (yy0 & 0b01) {
                                    continue;
                                }
                            }
                            _ => (),
                        }

                        rmv.push((
                            pos0..pos1,
                            [rid0, pos0, pos1, rid1, ppos0, ppos1, orientation, group],
                        ));
                        //DEBUG
                        //println!("{} {} {} {} {} {} {} {}", rid0, pos0, pos1, rid1, ppos0, ppos1, orientation, group);
                    }
                }
                rmv
            })
            .collect_into_vec(&mut range_matches);

        for rmv in range_matches {
            itvl_vec.extend(rmv);
        }

        let itvl: IntervalTree<u32, MapIntervalRecord> = itvl_vec.iter().cloned().collect();
        all_itvl.insert(rid0 as u32, itvl);
    }
    all_itvl
}

pub fn build_shmmer_map_from_query_results(mqr: &Vec<MapIntervalRecord>) -> MapIntervals {
    let mut id2itvls = FxHashMap::<u32, Vec<(Range<u32>, MapIntervalRecord)>>::default();

    mqr.iter().for_each(|&r| {
        id2itvls
            .entry(r[0])
            .or_insert_with(|| vec![])
            .push((r[1]..r[2], r));
    });

    let mut all_itvl = MapIntervals::default();
    id2itvls.iter().for_each(|(sid, recs)| {
        let itvl: IntervalTree<u32, MapIntervalRecord> = recs.iter().cloned().collect();
        all_itvl.insert(*sid, itvl);
    });

    all_itvl
}

pub fn map_interval_query(
    ivtl: &MapIntervals,
    sid: u32,
    bgn: u32,
    end: u32,
) -> Vec<MapIntervalRecord> {
    let mut q_res: Vec<_> = ivtl
        .get(&sid)
        .unwrap()
        .query(bgn..end)
        .map(|x| x.value)
        .collect();
    q_res.sort();
    let mut mq = Vec::<MapIntervalRecord>::new();
    for e in q_res {
        mq.push(e);
        // DBEUG
        // println!("q: {} {} {} r:{} {} {} {} {} {} {}", sid, bgn, end, e[0], e[1], e[2], e[3], e[4], e[5], e[6]);
    }
    mq
}

pub fn find_match_chain(matches: &Vec<MapIntervalRecord>) -> Vec<MapIntervalRecord> {
    let mut seqpair_count = FxHashMap::<(u32, u32), u32>::default();
    let mut out = Vec::<MapIntervalRecord>::new();
    matches.iter().for_each(|v| {
        let sid0 = v[0];
        let sid1 = v[3];
        *seqpair_count.entry((sid0, sid1)).or_insert(0) += 1;
    });

    for ((sid0, sid1), count) in seqpair_count.iter() {
        if *count < 4 {
            continue;
        }
        let mut aln_g = DiGraphMap::<u32, u32>::new();
        let filtered_matches = matches
            .iter()
            .filter(|x| x[0] == *sid0 && x[3] == *sid1)
            .collect::<Vec<&MapIntervalRecord>>();
        for i in 0..filtered_matches.len() {
            let mut j = i + 1;
            loop {
                if j >= filtered_matches.len() || j > i + 25 {
                    break;
                }
                let r1 = filtered_matches[i];
                let r2 = filtered_matches[j];
                if r1[1] == r2[1] {
                    j += 1;
                    continue;
                }
                let d1: i64;
                let d2: i64;
                if r1[6] == 1 {
                    d1 = r1[1] as i64 + r1[4] as i64;
                } else {
                    d1 = r1[1] as i64 - r1[4] as i64;
                }
                if r2[6] == 1 {
                    d2 = r2[1] as i64 + r2[4] as i64;
                } else {
                    d2 = r2[1] as i64 - r2[4] as i64;
                }
                let d: f32 =
                    ((d1 - d2).abs() as f32) / ((r2[1] as i64 - r1[1] as i64).abs() as f32);
                if d < 0.2 && r1[6] == r2[6] && r2[1] - r1[1] < 100000 {
                    aln_g.add_edge(i as u32, j as u32, 1);
                }
                j += 1;
            }
        }
        let mut vertex_sets = UnionFind::new(filtered_matches.len());
        for e in aln_g.all_edges() {
            vertex_sets.union(e.0, e.1);
        }

        let labels = vertex_sets.into_labeling();
        for i in 0..filtered_matches.len() {
            if !aln_g.contains_node(i as u32) {
                continue;
            }
            let label = labels[i];
            //println!("{} {} {:?} {} {}", sid0, sid1, filtered_matches[i], i, label);
            let mut out_r: [u32; 8] = filtered_matches[i].clone();
            out_r[7] = label;
            out.push(out_r);
        }
    }
    out
}

pub struct ExtDeltas {
    pub sid0: u32,
    pub pos0: u32,
    pub sid1: u32,
    pub pos1: u32,
    pub dk: i32,
    pub strand1: u8,
    pub base: char,
}

pub fn get_deltas(
    subseq0: &Vec<u8>,
    subseq1: &Vec<u8>,
    sid0: u32,
    sid1: u32,
    offset0: u32,
    offset1: u32,
    strand: u8,
) -> Vec<ExtDeltas> {
    let mut out_vec = Vec::<ExtDeltas>::new();
    if subseq0.len() < 32 || subseq1.len() < 32 {
        return out_vec;
    }
    //let hs0 = super::shmmrutils::get_hpc_seq(&subseq0).s;
    //let hs1 = super::shmmrutils::get_hpc_seq(&subseq1).s;
    let hs0 = subseq0;
    let hs1 = subseq1;

    if let Some(ovlp) = super::shmmrutils::match_reads(hs0, hs1, true, 0.02, 32, 24) {
        if let Some(mut deltas) = ovlp.deltas {
            if deltas.len() > 0 {
                deltas.reverse();
                for dpt in deltas {
                    let base = if dpt.dk > 0 {
                        '-'
                    } else {
                        subseq1[dpt.y as usize - 1] as char
                    };
                    out_vec.push(ExtDeltas {
                        sid0: sid0,
                        pos0: offset0 + dpt.x,
                        sid1: sid1,
                        pos1: offset1 + dpt.y,
                        strand1: strand,
                        dk: dpt.dk,
                        base: base,
                    });
                }
            }
        }
    }
    out_vec
}

fn generate_deltas(seqdb0: &SeqDB, seqdb1: &SeqDB, m: MapIntervalRecord) -> usize {
    let mut rev = false;
    if m[6] == 1 {
        //strand
        rev = true;
    }
    let sid0 = m[0] as usize;
    let b0 = m[1] as usize;
    let e0 = m[2] as usize;
    let sid1 = m[3] as usize;
    let b1 = m[4] as usize;
    let e1 = m[5] as usize;
    let subseq0: Vec<u8> = seqdb0.seqs[sid0 as usize][b0..e0]
        .iter()
        .map(|c| c & 0x0F)
        .collect();

    let subseq1: Vec<u8>;
    let strand: u8;

    match rev {
        true => {
            let bb1 = b1 - (56 - 2);
            let ee1 = e1 - (56 - 2);
            let slen = seqdb1.seqs[sid1].len();
            let b = slen - ee1;
            let e = slen - bb1;

            subseq1 = seqdb1.seqs[sid1 as usize][b..e]
                .iter()
                .map(|c| (c >> 4) & 0x0F)
                .collect();
            strand = 1;
            /*
            println!("DEBUB s0 {}", String::from_utf8(subseq0.clone()).unwrap());
            println!("DEBUB s1 {}", String::from_utf8(subseq1.clone()).unwrap());
            */
        }
        false => {
            subseq1 = seqdb1.seqs[sid1 as usize][b1..e1]
                .iter()
                .map(|c| c & 0x0F)
                .collect();
            strand = 0;
        }
    }
    let subseq0 = get_hpc_seq(&subseq0).s;
    let subseq1 = get_hpc_seq(&subseq1).s;

    let v = get_deltas(
        &subseq0,
        &subseq1,
        sid0 as u32,
        sid1 as u32,
        m[1],
        m[4],
        strand,
    );

    //let sname0 = &seqdb0.id2seqname[&(sid0 as u32)];
    //let sname1 = &seqdb1.id2seqname[&(sid1 as u32)];
    /*
    v.iter().for_each(|e| {
        println!(
            "D {} {} {} {} {} {} {} {} {}",
            sname0, e.sid0, e.pos0, sname1, e.sid1, e.pos1, e.strand1, e.dk, e.base
        );
    });
    */

    v.len()
}

pub fn linearize_hits(mut v: Vec<[u32; 8]>) -> Vec<[u32; 8]> {
    let mut vv = Vec::<[u32; 8]>::new();

    v.sort_by(|a, b| (a[1], a[3], a[4]).partial_cmp(&(b[1], b[3], b[4])).unwrap());
    let mut strand_count = FxHashMap::<u32, [u32; 2]>::default();
    v.iter().for_each(|x| {
        let sid = x[3];
        let strand = x[6];
        let sc = strand_count.entry(sid).or_insert([0_u32; 2]);
        if strand == 0 {
            sc[0] += 1;
        } else {
            sc[1] += 1;
        }
    });

    // identify the major alignment strand between two contigs
    let mut strands = FxHashMap::<u32, u32>::default();
    strand_count.iter().for_each(|(k, v)| {
        if v[0] < v[1] {
            strands.entry(*k).or_insert(1);
        } else {
            strands.entry(*k).or_insert(0);
        }
    });

    let mut last_postions = FxHashMap::<u32, [u32; 2]>::default();
    v.iter().for_each(|x| {
        let sid = x[3];
        if !last_postions.contains_key(&sid) {
            vv.push(x.clone());
            let strand = *strands.get(&sid).unwrap();
            if x[6] != strand {
                return;
            }
            if strand == 0 {
                last_postions.insert(sid, [x[2], x[5]]);
            } else {
                last_postions.insert(sid, [x[2], x[4]]);
            }
        } else {
            let last_p = last_postions.get(&sid).unwrap();
            let strand = *strands.get(&sid).unwrap();
            if x[6] != strand {
                return;
            }
            if x[2] <= last_p[0] {
                return;
            }
            if strand == 0 {
                if x[5] <= last_p[1] {
                    return;
                } else {
                    vv.push(x.clone());
                    last_postions.insert(sid, [x[2], x[5]]);
                }
            } else {
                if x[4] >= last_p[1] {
                    return;
                } else {
                    vv.push(x.clone());
                    last_postions.insert(sid, [x[2], x[4]]);
                }
            }
        }
    });
    vv
}

fn filter_chain_group(chain_groups: FxHashMap<(u32, u32, u32), Vec<[u32; 8]>>) -> Vec<[u32; 8]> {
    let mut chain_list = chain_groups
        .iter()
        .collect::<Vec<(&(u32, u32, u32), &Vec<[u32; 8]>)>>();

    chain_list.sort_by_key(|v| -(v.1.len() as i64));

    let mut r_intervals = FxHashSet::<(u32, u32, u32)>::default();
    let mut t_intervals = FxHashSet::<(u32, u32, u32)>::default();
    let mut aln_itvls = Vec::<[u32; 8]>::new();
    for (_k, v) in chain_list.iter() {
        //DEBUG
        //println!("# {:?} {:?}", _k, v.len());
        if v.len() < 4 {
            continue;
        }
        let mut aln_itvls0 = Vec::<[u32; 8]>::new();
        for w in v.iter() {
            let r_ivtl = (w[3], w[1], w[2]); // tagged with target id for both
            let t_ivtl = (w[3], w[4], w[5]);
            if !r_intervals.contains(&r_ivtl) && !t_intervals.contains(&t_ivtl) {
                aln_itvls0.push(*w);
                r_intervals.insert(r_ivtl);
                t_intervals.insert(t_ivtl);
            }
        }
        if aln_itvls0.len() > 0 && aln_itvls0.len() as f32 > 0.5 * v.len() as f32 {
            //ad hoc rule to avoid repeat
            aln_itvls.extend(aln_itvls0);
        }
    }
    aln_itvls
}

pub fn dedup_target_seqs(
    ref_fasta_file: &String,
    target_fasta_file: &String,
    output_file: &String,
    w: u32,
    k: u32,
    r: u32,
) -> Result<(), std::io::Error> {
    let mut sdb0 = SeqDB::new(ref_fasta_file.clone());
    let mut sdb1 = SeqDB::new(target_fasta_file.clone());

    sdb0.load_sequences()?;
    sdb1.load_sequences()?;
    sdb0.build_shmmrs(w, k, r);
    sdb1.build_shmmrs(w, k, r);

    let shmmrmap = generate_shmmr_map(&sdb0.shmmrs, &sdb1.shmmrs, 32);
    let mut out_file = BufWriter::new(File::create(format!("{}", output_file))?);

    let map_interval_list = (0..sdb0.seqs.len())
        .into_par_iter()
        .map(|sid| {
            let seqname = &sdb0.id2seqname[&(sid as u32)];
            let seqlen = sdb0.seqlen[seqname];
            let mq = map_interval_query(&shmmrmap, sid as u32, 0, seqlen as u32);
            let mq_filtered = find_match_chain(&mq);
            let mut chain_group = FxHashMap::<(u32, u32, u32), Vec<[u32; 8]>>::default();
            mq_filtered.iter().for_each(|r| {
                let (sid0, sid1, chain_label) = (r[0], r[3], r[7]);
                let key = (sid0, sid1, chain_label);
                //DEBUG
                /*
                println!(
                    "{} {} {} {} {} {} {} {} {}",
                    r[0],
                    r[1],
                    r[2],
                    r[3],
                    r[4],
                    r[5],
                    r[6],
                    r[7],
                    sdb0.id2seqname.get(&r[3]).unwrap()
                );
                */
                //println!("# {} {} {}", sid0, sid1, chain_label );
                chain_group.entry(key).or_insert_with(|| vec![]).push(*r);
            });

            let new_chains = filter_chain_group(chain_group);
            /*
            let new_chains = chain_group
                .values()
                .into_iter()
                .filter(|chain| chain.len() >= 0)
                .map(|chain| chain.clone())
                .flatten()
                .collect::<Vec<[u32; 8]>>();
                */
            (sid, seqlen, new_chains)
        })
        .filter(|(_sid, _seqlen, new_chains)| new_chains.len() > 0)
        .map(|(sid, seqlen, new_chains)| {
            let shmmrmap = build_shmmer_map_from_query_results(&new_chains);
            let mq = map_interval_query(&shmmrmap, sid as u32, 0, seqlen as u32);

            //TODO: we need to build a new chain group and for overlapped chains, we only use the longest one

            let v = mq
                .into_iter()
                .map(|m| {
                    let dcount = generate_deltas(&sdb0, &sdb1, m);
                    let mut out = [0_u32; 8];
                    out[0..7].copy_from_slice(&m[0..7]);
                    out[7] = dcount as u32;
                    out
                })
                .collect::<Vec<[u32; 8]>>();

            linearize_hits(v)
        })
        .flatten()
        .collect::<Vec<[u32; 8]>>();

    let mut map_intervals = FxHashMap::<[u32; 3], Vec<[u32; 8]>>::default();
    map_interval_list.iter().for_each(|r| {
        let key = [r[0], r[3], r[6]];
        map_intervals.entry(key).or_insert_with(|| vec![]).push(*r);
    });

    let mut dup_seq_ids = FxHashSet::<u32>::default();
    map_intervals.iter().for_each(|(k, vv)| {
        let mut aligned_base_count = 0_u32;
        let mut delta_count = 0_u32;
        let [sid0, sid1, _strand] = *k;
        vv.iter().for_each(|v| {
            aligned_base_count += v[5] - v[4];
            delta_count += v[7];
        });
        let seqlen0 = sdb0.seqs[sid0 as usize].len();
        let seqlen1 = sdb1.seqs[sid1 as usize].len();
        //println!("C {} {} {} {} {} {}", k[0], k[1], k[2], aligned_base_count, delta_count, sdb1.seqs[k[1] as usize].len());
        if delta_count == 0 && (seqlen1 as u32) - aligned_base_count <= 2000 {
            if seqlen0 == seqlen1 {
                // if there are two identifcal sequecnes, we only keep one
                if sid1 > sid0 {
                    dup_seq_ids.insert(sid1);
                }
            } else {
                dup_seq_ids.insert(sid1);
            }
        };
    });

    (0..sdb1.seqs.len()).try_for_each(|sid1| -> Result<(), std::io::Error> {
        let sid1 = sid1 as u32;
        let base_map = &[b'A', b'C', b'G', b'T', b'a', b'c', b'g', b't'];
        if !dup_seq_ids.contains(&sid1) {
            let seq_name = sdb1.id2seqname.get(&sid1).unwrap();
            if sdb1.masks[sid1 as usize].len() != 1 {
                return Ok(());
            }
            writeln!(out_file, ">{}", seq_name)?;
            let mut seq_tmp: Vec<u8> = sdb1.seqs[sid1 as usize]
                .iter()
                .map(|b| (b & 0b0011) + 4)
                .collect();
            sdb1.masks[sid1 as usize].iter().for_each(|&(b, e)| {
                for i in b..e {
                    seq_tmp[i] -= 4;
                }
            });
            let seq: Vec<u8> = seq_tmp
                .iter()
                .filter(|&b| *b < 4)
                .map(|&b| base_map[b as usize])
                .collect();
            writeln!(out_file, "{}", String::from_utf8_lossy(&seq))?;
        }
        Ok(())
    })?;

    Ok(())
}
