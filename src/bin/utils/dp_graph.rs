// Peregrine Assembler and SHIMMER Genome Assembly Toolkit 
// 2019, 2020, 2021- (c) by Jason, Chen-Shan, Chin
//
// This Source Code Form is subject to the terms of the 
// Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
//
// You should have received a copy of the license along with this
// work. If not, see <http://creativecommons.org/licenses/by-nc-sa/4.0/>.

#![allow(dead_code)]

///
/// This code convert  the overlaps into assembly graph. Like most genome 
/// assembler, there are many huristici rules guided by some general 
/// principles and intuitions.
/// 

use glob::glob;
use petgraph::graphmap::DiGraphMap;
use petgraph::visit::{Bfs, DfsPostOrder};
use petgraph::Direction::{Incoming, Outgoing};
use rustc_hash::{FxHashMap, FxHashSet};

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::thread;

use super::graph_analysis::*;
use super::log_resource;
use super::{getrusage, MaybeUninit, RUSAGE_SELF};
use std::io::prelude::*;

fn build_read_ovlp_data<P>(filename: P) -> OverlapMap
where
    P: AsRef<Path>,
{
    // 
    // parsing the overlapper output converting into internal data structure `OverlapMap`
    // 

    let mut rid2ovlp = OverlapMap::default();
    let mut buffer = String::new();

    let file = File::open(filename);
    let _err: Result<usize, io::Error> = file.unwrap().read_to_string(&mut buffer);
    for line in buffer.split("\n") {
        let mut v: Vec<&str> = Vec::<&str>::with_capacity(24); // we need pre-allocate some space for performance
        line.split(' ').for_each(|c| v.push(c));
        match v[0] {
            "O" => {
                let ovlp = Overlap::build_from(v);

                if ovlp.dist >= 12 && ovlp.max_dist_c == 0 {
                    // if overapped reads have big differences in homoployer regions (dist large but dist_c smalle)
                    // we might need to improve this in lower coverage cases 
                    continue;
                }

                if ovlp.dist_c > 2 {
                    // likely to be true differences
                    continue;
                }

                rid2ovlp
                    .entry(ovlp.rid0)
                    .or_insert_with(|| vec![])
                    .push(ovlp);
            }
            _ => (),
        }
    }
    rid2ovlp
}

fn is_dead_ended(v: &Vec<Overlap>, s: &FxHashSet<u32>) -> bool {
    // check if the read only have one end that is overlapped with others
    let mut right_ext = false;
    let mut left_ext = false;
    for vv in v.iter() {
        if s.contains(&vv.rid1) {
            continue;
        }

        if vv.d_right > 0 {
            right_ext = true;
        }
        if vv.d_left < 0 {
            left_ext = true;
        }
    }
    if right_ext && left_ext {
        false
    } else {
        true
    }
}

fn is_chimer(v: &Vec<Overlap>) -> bool {
    // huristic method to determine if a read may be a chimer by a quick coverage over the read
    // analysis 
    let mut left_most: i32 = i32::MAX;
    let mut right_most: i32 = i32::MIN;
    let mut rlen: u32 = 0;
    for vv in v.iter() {
        //if (vv.max_dist_c > 4 && vv.dist_c > 2) || (vv.max_dist_c <= 4 && vv.dist_c > 1) {
        //    continue;
        //}
        if vv.dist_c > 2 {
            continue;
        }
        if rlen == 0 {
            rlen = vv.len0;
        }
        let d_left = vv.d_left;
        let d_right = vv.d_right;
        if d_right > 0 && d_left < left_most {
            left_most = d_left;
        }
        if d_left < 0 && d_right > right_most {
            right_most = d_right;
        }
    }
    if rlen == 0 || (rlen as i32) + right_most < left_most {
        true
    } else {
        false
    }
}

fn get_rpair2ovlps(
    rid2ovlp_all: &FxHashMap<u32, Vec<Overlap>>,
    contained: &FxHashSet<u32>,
    chimers: &FxHashSet<u32>,
    low_q: &FxHashSet<u32>,
    bestn: usize,
) -> FxHashMap<ReadPair, Overlap> {

    // filter overlaps to find all dovetail overlaps 
    // it needs the informaiton of contained / chimers / low_q reads etc, to 
    // avoid missing overlaps masked by those.
    //
    // we also reduce the overlaps used to only the best reads to reduce computation 
    // complexity.

    let mut rpair2overlap = FxHashMap::<ReadPair, Overlap>::default();
    let mut all_rids = Vec::<u32>::new();
    for r in rid2ovlp_all.keys() {
        all_rids.push(*r);
    }
    all_rids.sort();
    for r in all_rids {
        let v = rid2ovlp_all.get(&r).unwrap();
        let mut left_candidates = Vec::<Overlap>::with_capacity(64);
        let mut right_candidates = Vec::<Overlap>::with_capacity(64);

        if contained.contains(&r) {
            // println!("T0 {} {}", r, v.len());
            continue;
        }

        if chimers.contains(&r) {
            // println!("T0 {} {}", r, v.len());
            continue;
        }

        if low_q.contains(&r) {
            continue;
        }

        for vv in v.iter() {
            if contained.contains(&vv.rid1) {
                continue;
            }
            if chimers.contains(&vv.rid1) {
                continue;
            }
            if low_q.contains(&vv.rid1) {
                continue;
            }
            // let _ovlp_len = vv.end0 - vv.bgn0;
            if vv.d_right > 0 && vv.d_left > 0 {
                right_candidates.push(*vv);
            }
            if vv.d_right < 0 && vv.d_left < 0 {
                left_candidates.push(*vv);
            }
        }
        /*
        log::debug!(
            "get_rpair2ovlps: T {} {} {} {}",
            r,
            v.len(),
            left_candidates.len(),
            right_candidates.len()
        );
        */
        left_candidates.sort_by(|a, b| {
            let la = a.end0 - a.bgn0;
            let lb = b.end0 - b.bgn0;
            lb.cmp(&la)
        });

        let mut found = false;

        for i in 0..bestn {
            if let Some(vv) = left_candidates.get(i) {
                if vv.dist == 0 {
                    // println!("O {}", vv.format());
                    let rp = ReadPair {
                        rid0: vv.rid1,
                        strand0: vv.strand1,
                        rid1: vv.rid0,
                        strand1: 0,
                    };
                    rpair2overlap.insert(rp, vv.swap_rp());
                    let rp = rp.reverse();
                    rpair2overlap.insert(rp, *vv);
                    found = true;
                }
            }
        }

        if !found {
            for i in 0..bestn {
                if let Some(vv) = left_candidates.get(i) {
                    if vv.dist_c == 0 {
                        // println!("O {}", vv.format());
                        let rp = ReadPair {
                            rid0: vv.rid1,
                            strand0: vv.strand1,
                            rid1: vv.rid0,
                            strand1: 0,
                        };
                        rpair2overlap.insert(rp, vv.swap_rp());
                        let rp = rp.reverse();
                        rpair2overlap.insert(rp, *vv);
                        found = true;
                    }
                }
            }
        }

        if !found {
            for i in 0..bestn {
                if let Some(vv) = left_candidates.get(i) {
                    if vv.dist_c < 5 {
                        // println!("O {}", vv.format());
                        let rp = ReadPair {
                            rid0: vv.rid1,
                            strand0: vv.strand1,
                            rid1: vv.rid0,
                            strand1: 0,
                        };
                        rpair2overlap.insert(rp, vv.swap_rp());
                        let rp = rp.reverse();
                        rpair2overlap.insert(rp, *vv);
                    }
                }
            }
        }

        right_candidates.sort_by(|a, b| {
            let la = a.end0 - a.bgn0;
            let lb = b.end0 - b.bgn0;
            lb.cmp(&la)
        });

        // huristic to pick the best reads as the "right" overlaps
        // it falls back to slightly worse reads if necessary

        found = false;

        for i in 0..bestn {
            if let Some(vv) = right_candidates.get(i) {
                if vv.dist == 0 {
                    // println!("O {}", vv.format());
                    let rp = ReadPair {
                        rid0: vv.rid0,
                        strand0: 0,
                        rid1: vv.rid1,
                        strand1: vv.strand1,
                    };
                    rpair2overlap.insert(rp, *vv);
                    let rp = rp.reverse();
                    rpair2overlap.insert(rp, vv.swap_rp());
                    found = true;
                }
            }
        }

        if !found {
            for i in 0..bestn {
                if let Some(vv) = right_candidates.get(i) {
                    if vv.dist_c == 0 { // using dist_c than dist for hp compressed match
                        // println!("O {}", vv.format());
                        let rp = ReadPair {
                            rid0: vv.rid0,
                            strand0: 0,
                            rid1: vv.rid1,
                            strand1: vv.strand1,
                        };
                        rpair2overlap.insert(rp, *vv);
                        let rp = rp.reverse();
                        rpair2overlap.insert(rp, vv.swap_rp());
                        found = true;
                    }
                }
            }
        }

        if !found {
            for i in 0..bestn {
                if let Some(vv) = right_candidates.get(i) {
                    // println!("O {}", vv.format());
                    if vv.dist_c < 5 {
                        let rp = ReadPair {
                            rid0: vv.rid0,
                            strand0: 0,
                            rid1: vv.rid1,
                            strand1: vv.strand1,
                        };
                        rpair2overlap.insert(rp, *vv);
                        let rp = rp.reverse();
                        rpair2overlap.insert(rp, vv.swap_rp());
                    }
                }
            }
        }
    }
    rpair2overlap
}

fn get_g0(rpair2overlap: &FxHashMap<ReadPair, Overlap>) -> DiGraphMap<(u32, u8), u32> {

    // get `g0` the raw assembly overlap graph

    let mut g0 = DiGraphMap::<(u32, u8), u32>::new();
    for rp in rpair2overlap.keys() {
        let rid0 = rp.rid0;
        let strand0 = rp.strand0;
        let rid1 = rp.rid1;
        let strand1 = rp.strand1;
        let ovlp = rpair2overlap.get(rp).unwrap();
        let ovlp_len = ovlp.end0 - ovlp.bgn0;
        g0.add_edge((rid0, strand0), (rid1, strand1), ovlp_len);
        //println!("E {} {}", rp.to_str(), ovlp_len)
    }
    g0
}

fn patch_ends(
    mut g0: &mut DiGraphMap<(u32, u8), u32>,
    _contained: &FxHashSet<u32>,
    rid2ovlp_all: &FxHashMap<u32, Vec<Overlap>>,
    rpair2overlap: &mut FxHashMap<ReadPair, Overlap>,
) -> () {

    // rescue the connection between two regions that are connected by a tangled blob

    let mut rdata = unsafe { MaybeUninit::uninit().assume_init() };
    let _res = unsafe { getrusage(RUSAGE_SELF, &mut rdata) };

    let mut bgn_nodes = FxHashSet::<(u32, u8)>::default();
    let mut end_nodes = FxHashSet::<(u32, u8)>::default();

    log_resource("BGN: patch_ends, stage1", &mut rdata);
    for v in g0.nodes() {
        let in_deg = g0.neighbors_directed(v, Incoming).count();
        let out_deg = g0.neighbors_directed(v, Outgoing).count();
        if in_deg == 0 && out_deg != 0 {
            bgn_nodes.insert(v);
            end_nodes.insert((v.0, 1 - v.1));
        }
    }
    log::info!("bgn_nodes len: {}", bgn_nodes.len());
    log_resource("END: patch_ends, stage1", &mut rdata);

    log_resource("BGN: patch_ends, stage2", &mut rdata);
    // handle some false containment cases
    let mut end_nodes_vec = end_nodes.into_iter().collect::<Vec<(u32, u8)>>();
    end_nodes_vec.sort();
    for r in end_nodes_vec.iter() {
        if !rid2ovlp_all.contains_key(&r.0) {
            continue;
        }

        let v = rid2ovlp_all.get(&r.0).unwrap();
        for vv in v.iter() {
            let rid0 = vv.rid0;
            let strand0 = 0;
            let rid1 = vv.rid1;
            let strand1 = vv.strand1;
            let ovlp_len = vv.end0 - vv.bgn0;

            if vv.dist >= 12 && vv.max_dist_c == 0 {
                continue;
            }

            if vv.d_right > 0 && vv.d_left > 0 {
                g0.add_edge((rid0, strand0), (rid1, strand1), ovlp_len);
                g0.add_edge((rid1, 1 - strand1), (rid0, 1 - strand0), ovlp_len);
                let rp = ReadPair {
                    rid0: vv.rid0,
                    strand0: 0,
                    rid1: vv.rid1,
                    strand1: vv.strand1,
                };
                rpair2overlap.insert(rp, *vv);
                let rp = rp.reverse();
                rpair2overlap.insert(rp, vv.swap_rp());
            }
            if vv.d_right < 0 && vv.d_left < 0 {
                g0.add_edge((rid1, strand1), (rid0, strand0), ovlp_len);
                g0.add_edge((rid0, 1 - strand0), (rid1, 1 - strand1), ovlp_len);
                let rp = ReadPair {
                    rid0: vv.rid1,
                    strand0: vv.strand1,
                    rid1: vv.rid0,
                    strand1: 0,
                };
                rpair2overlap.insert(rp, vv.swap_rp());
                let rp = rp.reverse();
                rpair2overlap.insert(rp, *vv);
            }
        }
    }
    log_resource("END: patch_ends, stage2", &mut rdata);

    log_resource("BGN: patch_ends, stage3", &mut rdata);
    transitive_reduction(&mut g0);
    log_resource("END: patch_ends, stage3", &mut rdata);
    
    log_resource("BGN: patch_ends, stage5", &mut rdata);
    transitive_reduction(&mut g0);
    log_resource("END: patch_ends, stage5", &mut rdata);
}

fn is_branch(g: &DiGraphMap<(u32, u8), u32>, w: &(u32, u8), path_reads: &FxHashSet<u32>) -> bool {

    // Given a set of reads in a path and a node, check if there is ture branch at the node in the paths
    // If all out edge eventually merge, then we see a bubble rather than a true branch.

    let mut out_nodes = Vec::<(u32, u8)>::new();
    let mut branched = false;
    if g.neighbors_directed(*w, Outgoing).count() > 1 {
        for u in g.neighbors_directed(*w, Outgoing) {
            if path_reads.contains(&u.0) {
                //handle invert repeat cases
                continue;
            }
            out_nodes.push(u);
        }

        let out_count = out_nodes.len();
        if out_count == 0 {
            return branched;
        }

        for u in out_nodes {
            let mut bfs = Bfs::new(g, u);
            let mut alt_path_len = 0_u32;
            loop {
                if let Some(n) = bfs.next(g) {
                    alt_path_len += 1;
                    if path_reads.contains(&n.0) {
                        break;
                    }
                    if alt_path_len > 1024 {
                        branched = true;
                        break;
                    }
                } else {
                    break;
                }
            }
            //log::debug!("outnode: {:?} {:?} {} {:?}", w, u, alt_path_len, branched);
            if branched {
                break;
            }
        }
    }
    branched
}

fn trim_path(
    g: &DiGraphMap<(u32, u8), u32>,
    path_edge_list: &Vec<((u32, u8), (u32, u8))>,
) -> Vec<((u32, u8), (u32, u8))> {

    // 
    // given a list of edges, generate an unbranch path from the beginning
    //

    let mut out_edge_list = Vec::<((u32, u8), (u32, u8))>::new();
    let mut path_reads = FxHashSet::<u32>::default();
    for &(v, w) in path_edge_list.iter() {
        if path_reads.len() == 0 {
            path_reads.insert(v.0);
        }
        path_reads.insert(w.0);
    }

    for &(v, w) in path_edge_list.iter() {
        out_edge_list.push((v, w));
        if is_branch(g, &w, &path_reads) {
            break;
        }

        if is_branch(g, &(w.0, 1 - w.1), &path_reads) {
            break;
        }
    }
    out_edge_list
}

fn get_path_from_seed(
    g1: &DiGraphMap<(u32, u8), u32>,
    seed: &(u32, u8),
) -> (Vec<((u32, u8), (u32, u8))>, i32) {
  
    // find the best path from a seed


    let v = *seed;
    let mut path_weight = FxHashMap::<(u32, u8), i32>::default();

    let mut max_pre = FxHashMap::<Option<(u32, u8)>, Option<(u32, u8)>>::default();
    path_weight.insert((v.0, v.1), 0);

    let mut best_weight_end = Option::<(i32, u32, u8)>::None;
    // count the number of reads used in the path
    let mut best_score = 0_i32;

    let mut post_order = DfsPostOrder::new(&g1, v);
    let mut rev_post_order = Vec::<(u32, u8)>::new();
    loop {
        if let Some(n) = post_order.next(&g1) {
            rev_post_order.push(n);
            log::debug!("layout debug: n {:?} {:?}", seed, n);
        } else {
            break;
        }
    }
    rev_post_order.reverse();

    let mut read_count = FxHashMap::<u32, u32>::default();
    for n in rev_post_order {
        let mut m = 0_i32;
        let mut p = Option::<(u32, u8)>::None;
        for w in g1.neighbors_directed(n, Incoming) {
            if let Some(m0) = path_weight.get(&w) {
                if *m0 > m {
                    m = *m0;
                    p = Some(w);
                }
            }
        }

        if !read_count.contains_key(&n.0) {
            path_weight.insert(n, m + 1);
            if m + 1 > best_score {
                best_weight_end = Some((m + 1, n.0, n.1));
                best_score = m + 1;
            }
        } else {
            path_weight.insert(n, m - 1); // if there is invert repeat, reduce the weight
        }

        max_pre.insert(Some(n), p);

        let entry = read_count.entry(n.0).or_insert(0);
        *entry += 1;
    }

    let mut path_edge_list = Vec::<((u32, u8), (u32, u8))>::new();

    if let Some((_, w0, w1)) = best_weight_end {
        let mut w = Some((w0, w1));
        let mut used_edges = FxHashSet::<((u32, u8), (u32, u8))>::default();
        loop {
            match max_pre.get(&w).unwrap() {
                Some(v) => {
                    let e = (*v, w.unwrap());
                    path_edge_list.push(e);
                    if used_edges.contains(&e) {
                        break;
                    } else {
                        used_edges.insert(e);
                    }
                    w = Some(*v);
                }
                None => break,
            }
        }
        path_edge_list.reverse();
    }

    let mut path_reads = FxHashSet::<u32>::default();
    for &(v, w) in path_edge_list.iter() {
        if path_reads.len() == 0 {
            path_reads.insert(v.0);
        }
        path_reads.insert(w.0);
    }

    // we don't update the best_score here after trimming so we won't create long inverted repeated contig
    let out_edge_list = trim_path(g1, &path_edge_list);
    (out_edge_list, best_score)
}

fn write_paths_from_seeds(
    g1: &mut DiGraphMap<(u32, u8), u32>,
    seeds: &Vec<(u32, u8, Vec<((u32, u8), (u32, u8))>, i32)>,
    used_edges: &mut FxHashSet<((u32, u8), (u32, u8))>,
    layout_file: &mut BufWriter<File>,
    rpair2overlap: &FxHashMap<ReadPair, Overlap>,
    ctg_id: &mut u32,
    ctg_tags: &mut FxHashMap<u32, &str>,
    read_to_ctg: &mut FxHashMap<ReadNode, Vec<(u32, u32)>>,
) -> () {

    // output the path from seeds in the layout file

    for v in seeds {
        let edge_list = &v.2;

        let mut ctg_tag = "P";
        if edge_list.len() < 3 {
            continue;
        }
        if edge_list.len() < 8 {
            // output some of these debris somewhere sometime leater
            // we might want to find the root cause of these short piece and properly proceess them later

            let mut attached = false;
            let b = edge_list.get(0).unwrap().0;
            //let e = edge_list.get(edge_list.len() - 1).unwrap().1;
            let mut bgn_contigs = FxHashSet::<u32>::default();
            if let Some(ns) = read_to_ctg.get(&b) {
                for n in ns {
                    if *ctg_tags.get(&n.0).unwrap() != "P" {
                        continue;
                    }
                    bgn_contigs.insert(n.0);
                }
            }
            if let Some(ns) = read_to_ctg.get(&(b.0, 1 - b.1)) {
                for n in ns {
                    if *ctg_tags.get(&n.0).unwrap() != "P" {
                        continue;
                    }
                    bgn_contigs.insert(n.0);
                }
            }

            for (_vv, ww) in edge_list {
                if let Some(ns) = read_to_ctg.get(ww) {
                    for n in ns {
                        if bgn_contigs.contains(&n.0) {
                            attached = true;
                            ctg_tag = "A";
                            break;
                        }
                    }
                }
                if attached {
                    ctg_tag = "A";
                    break;
                };
                if let Some(ns) = read_to_ctg.get(&(ww.0, 1 - ww.1)) {
                    for n in ns {
                        if bgn_contigs.contains(&n.0) {
                            attached = true;
                            ctg_tag = "A";
                            break;
                        }
                    }
                }
                if attached {
                    ctg_tag = "A";
                    break;
                };
            }

            if !attached {
                ctg_tag = "D"
            }
        }

        let (_, w) = edge_list[edge_list.len() - 1];
        let _res = writeln!(
            layout_file,
            "{} {} {}:{} {} {}:{} {}",
            ctg_tag,
            ctg_id,
            v.0,
            v.1,
            v.2.len(),
            w.0,
            w.1,
            edge_list.len()
        );

        let mut ctg_pos = 0_u32;
        for &(v, w) in edge_list {

            used_edges.insert((v, w));
            g1.remove_edge(v, w);
            g1.remove_edge((w.0, 1 - w.1), (v.0, 1 - v.1));

            let rp = ReadPair {
                rid0: v.0,
                strand0: v.1,
                rid1: w.0,
                strand1: w.1,
            };
            let mut ovlp = *rpair2overlap.get(&rp).unwrap();
            if v.1 == 1 {
                ovlp = ovlp.reverse_strand();
            }
            let _res = writeln!(
                layout_file,
                "E {} {} {} {} {} {} {} {} {}",
                ctg_id, v.0, v.1, w.0, w.1, ovlp.bgn0, ovlp.bgn1, ovlp.dist, ovlp.dist_c
            );
            let ovlp_len = ovlp.end0 - ovlp.bgn0;

            read_to_ctg
                .entry(v)
                .or_insert_with(|| vec![])
                .push((ctg_id.clone(), ctg_pos));

            ctg_pos += ovlp.len1 - ovlp_len;

            read_to_ctg
                .entry(w)
                .or_insert_with(|| vec![])
                .push((ctg_id.clone(), ctg_pos));
        }
        ctg_tags.insert(*ctg_id, ctg_tag);
        *ctg_id += 1;
    }

    let mut nodes_to_remove = Vec::<(u32, u8)>::new();
    g1.nodes()
        .into_iter()
        .filter(|&v| {
            g1.neighbors_directed(v, Incoming).count() == 0
                && g1.neighbors_directed(v, Outgoing).count() == 0
        })
        .for_each(|v| nodes_to_remove.push(v));
    nodes_to_remove.into_iter().for_each(|v| {
        g1.remove_node(v);
    });
}

fn generate_layout(
    g: &DiGraphMap<(u32, u8), u32>,
    rpair2overlap: &FxHashMap<ReadPair, Overlap>,
    mut read_to_ctg: &mut FxHashMap<ReadNode, Vec<(u32, u32)>>,
    output_prefix: &str,
) -> () {

    //
    // Interative generate the layour from the best seeds in the part of assembly without a circle.
    // After that, we resolve the sub graph with circles.
    //

    let layout_filename = format!("{}_layout.dat", output_prefix);
    let mut layout_file = BufWriter::new(File::create(layout_filename).unwrap());
    let mut ctg_id = 0_u32;

    let mut wg = g.clone();

    let mut ctg_tags = FxHashMap::<u32, &str>::default();
    loop {
        log::debug!("layout debug: in layout loop");
        let mut seeds = Vec::<(u32, u8, Vec<((u32, u8), (u32, u8))>, i32)>::new(); // a map from each component to the seed nodes

        for v in wg.nodes() {
            let in_degree = wg.neighbors_directed(v, Incoming).count();
            let out_degree = wg.neighbors_directed(v, Outgoing).count();

            if in_degree == 0 && out_degree != 0 {
                let (edge_list, score) = get_path_from_seed(&wg, &v);
                log::debug!("layout debug:  {:?} {:?}", v, edge_list.len());
                seeds.push((v.0, v.1, edge_list, score));
            }
        }
        log::debug!("layout debug: in layout loop 2");
        if seeds.len() == 0 {
            break;
        }

        let mut used_edges = FxHashSet::<((u32, u8), (u32, u8))>::default();
        seeds.sort_by(|a, b| (b.3.partial_cmp(&a.3)).unwrap());

        log::debug!("layout debug: in layout loop 3");
        let mut seeds2 = Vec::<(u32, u8, Vec<((u32, u8), (u32, u8))>, i32)>::new();

        let mut edges = FxHashSet::<((u32, u8), (u32, u8))>::default();
        for s in seeds {
            // println!("debug: {}:{} {} {}", s.0, s.1, s.2.len(), s.3);
            let mut overlapped = false;
            for e in s.2.iter() {
                if edges.contains(&e) {
                    overlapped = true;
                    break;
                }
            }

            if !overlapped {
                for &(v, w) in s.2.iter() {
                    edges.insert((v, w));
                    edges.insert(((w.0, 1 - w.1), (v.0, 1 - v.1)));
                }
                seeds2.push(s);
            }
        }
        
        write_paths_from_seeds(
            &mut wg,
            &seeds2,
            &mut used_edges,
            &mut layout_file,
            &rpair2overlap,
            &mut ctg_id,
            &mut ctg_tags,
            &mut read_to_ctg,
        );

        if used_edges.len() == 0 {
            break;
        }
    }

    // process the rest un-layouted nodes, likely in circles

    let mut seeds = Vec::<(u32, u8, Vec<((u32, u8), (u32, u8))>, i32)>::new(); // a map from each component to the seed nodes
    let mut reachable = FxHashSet::<(u32, u8)>::default();
    for v in wg.nodes() {
        if reachable.contains(&v) {
            continue;
        }

        let (edge_list, score) = get_path_from_seed(&wg, &v);
        //println!("test: {:?} {}", v, edge_list.len());
        reachable.insert(v);
        for (_v, w) in edge_list.iter() {
            reachable.insert(*w);
        }
        seeds.push((v.0, v.1, edge_list, score));
    }

    seeds.sort_by(|a, b| (b.3.partial_cmp(&a.3)).unwrap());

    let mut new_seeds = Vec::<(u32, u8, Vec<((u32, u8), (u32, u8))>, i32)>::new(); // a map from each component to the seed nodes
    let mut used_edges = FxHashSet::<((u32, u8), (u32, u8))>::default();
    for (v0, v1, edge_list, score) in seeds {
        let mut new_list = Vec::<((u32, u8), (u32, u8))>::new();
        for e in edge_list {
            let (v, w) = e;
            if !used_edges.contains(&e) {
                new_list.push(e);
                used_edges.insert(e);
                used_edges.insert(((w.0, 1 - w.1), (v.0, 1 - v.1)));
            } else {
                break;
            }
        }
        new_seeds.push((v0, v1, new_list, score));
    }

    let mut used_edges = FxHashSet::<((u32, u8), (u32, u8))>::default();
    write_paths_from_seeds(
        &mut wg,
        &new_seeds,
        &mut used_edges,
        &mut layout_file,
        &rpair2overlap,
        &mut ctg_id,
        &mut ctg_tags,
        &mut read_to_ctg,
    );
}

pub fn ovlp2layout_v2(prefix: &String, out_prefix: &String, bestn: usize) -> Result<(), io::Error> {
    
    //
    // The whole overlap to layout workflow  
    //

    let prefix = prefix.clone();
    let infile_pattern = [prefix, "*".to_string()].concat();

    let mut rdata = unsafe { MaybeUninit::uninit().assume_init() };
    let _res = unsafe { getrusage(RUSAGE_SELF, &mut rdata) };
    log_resource("BGN: load overlaps", &mut rdata);

    let mut children = Vec::new();
    let mut _chunk: u8 = 0;
    for entry in glob(&infile_pattern).expect("Failed to read glob pattern") {
        match entry {
            Ok(path) => {
                let child = thread::spawn(move || {
                    let rid2ovlp = build_read_ovlp_data(path);
                    rid2ovlp
                });
                children.push(child);
            }
            Err(e) => println!("{:?}", e),
        }
        _chunk += 1;
    }

    let mut rid2ovlp_all = OverlapMap::default();
    for child in children {
        let rid2ovlp_p = child.join().expect("oops! the child thread panicked");
        rid2ovlp_all.extend(rid2ovlp_p);
    }

    // filter out low quality reads
    let low_q = FxHashSet::<u32>::default();

    let mut dead_ends = FxHashSet::<u32>::default();
    let mut chimers = FxHashSet::<u32>::default();
    for i in 0..2 {
        for r in rid2ovlp_all.keys() {
            let v = rid2ovlp_all.get(r).unwrap();
            if is_dead_ended(v, &dead_ends) {
                dead_ends.insert(*r);
            }
            if i == 0 {
                if is_chimer(v) {
                    chimers.insert(*r);
                }
            }
        }
    }

    let mut contained = FxHashSet::<u32>::default();
    for r in rid2ovlp_all.keys() {
        let v = rid2ovlp_all.get(r).unwrap();
        for vv in v.iter() {

            if vv.dist_c > 2 {
                continue;
            }

            if dead_ends.contains(&vv.rid0) || dead_ends.contains(&vv.rid1) {
                continue;
            }
            if vv.d_left > 0 && vv.d_right < 0 {
                contained.insert(vv.rid1);
                // println!("CN1 {}", vv.format());
            }
            if vv.d_left < 0 && vv.d_right > 0 {
                contained.insert(vv.rid0);
                // println!("CN0 {}", vv.format());
            }
        }
    }
    log_resource("END: load overlaps", &mut rdata);

    // graph processing
    log_resource("BGN: build overlaps", &mut rdata);
    let mut rpair2overlap = get_rpair2ovlps(&rid2ovlp_all, &contained, &chimers, &low_q, bestn);
    log_resource("END: build overlaps", &mut rdata);

    log_resource("BGN: build g0", &mut rdata);
    let mut g0 = get_g0(&rpair2overlap);
    transitive_reduction(&mut g0);
    remove_simple_spur(&mut g0, 0);
    remove_single_bridge(&mut g0, 1000000);
    let gout_filename = format!("{}_g0.dat", out_prefix);
    let mut graph0_file = BufWriter::new(File::create(gout_filename).unwrap());

    for (v, w, len) in g0.all_edges() {
        let _res = writeln!(graph0_file, "G {} {} {} {} {}", v.0, v.1, w.0, w.1, len);
    }
    graph0_file.flush().expect("file write error");
    log_resource("END: build g0", &mut rdata);

    log_resource("BGN: patch ends", &mut rdata);
    patch_ends(&mut g0, &contained, &rid2ovlp_all, &mut rpair2overlap);
    log_resource("END: patch ends", &mut rdata);

    log_resource("BGN: utg0", &mut rdata);
    let paths = get_utg_paths(&g0);
    let (utg_g0, g_out) = utg_reduction(&paths, &g0);
    let utg_filename = format!("{}_utg0.dat", out_prefix);
    dump_utg_paths(&paths, &utg_g0, &utg_filename)?;
    log_resource("END: utg0", &mut rdata);

    log_resource("BGN: layout", &mut rdata);
    let mut read_to_ctg = FxHashMap::<ReadNode, Vec<(u32, u32)>>::default();
    generate_layout(&g_out, &rpair2overlap, &mut read_to_ctg, &out_prefix);

    let gout_filename = format!("{}_g.dat", out_prefix);
    let mut graph1_file = BufWriter::new(File::create(gout_filename).unwrap());
    
    for (v, w, len) in g_out.all_edges() {
        let _res = writeln!(graph1_file, "G {} {} {} {} {}", v.0, v.1, w.0, w.1, len);
    }
    log_resource("END: layout", &mut rdata);

    log_resource("BGN: utg1", &mut rdata);
    let paths = get_utg_paths(&g_out);
    let (utg_g, _) = utg_reduction(&paths, &g_out);
    let utg_filename = format!("{}_utg.dat", out_prefix);
    dump_utg_paths(&paths, &utg_g, &utg_filename)?;
    let gfa_filename = format!("{}_utg.gfa", out_prefix);
    dump_utg_gfa(&paths, &utg_g, &rpair2overlap, &read_to_ctg, &gfa_filename)?;
    log_resource("END: utg1", &mut rdata);

    Ok(())
}
