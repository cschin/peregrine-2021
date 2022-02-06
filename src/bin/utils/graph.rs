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
// old graph data structures and layout code for ovlp2layout_v1
//


use glob::glob;
use petgraph::graphmap::DiGraphMap;
use petgraph::visit::Bfs;
use petgraph::visit::Dfs;
use petgraph::Direction::{Incoming, Outgoing};
use rustc_hash::FxHashMap;
use rustc_hash::FxHashSet;

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::thread;

use std::io::prelude::*;

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

#[derive(Debug, Copy, Clone, Hash, Eq, PartialEq)]
struct ReadPair {
    rid0: u32,
    strand0: u8,
    rid1: u32,
    strand1: u8,
}

impl ReadPair {
    fn _to_str(&self) -> String {
        format!(
            "{} {} {} {}",
            self.rid0, self.strand0, self.rid1, self.strand1
        )
    }

    fn reverse(&self) -> ReadPair {
        ReadPair {
            rid0: self.rid1,
            strand0: 1 - self.strand1,
            rid1: self.rid0,
            strand1: 1 - self.strand0,
        }
    }
}

impl Overlap {
    fn _new() -> Self {
        Self {
            rid0: 0,
            rid1: 0,
            strand1: 0,
            len0: 0,
            len1: 0,
            d_left: 0,
            d_right: 0,
            bgn0: 0,
            end0: 0,
            bgn1: 0,
            end1: 0,
            dist: 0,
            idt: 0.0,
            dist_c: 0,
            max_dist_c: 0,
            idt_c: 0.0,
            flag: 0,
        }
    }

    fn build_from(v: Vec<&str>) -> Self {
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

    fn _format(&self) -> String {
        format!(
            "{} {} {} {} {} {} {} {} {} {} {} {} {:.2} {} {} {:.2} {}",
            self.rid0,
            self.rid1,
            self.strand1,
            self.len0,
            self.len1,
            self.d_left,
            self.d_right,
            self.bgn0,
            self.end0,
            self.bgn1,
            self.end1,
            self.dist,
            self.idt,
            self.dist_c,
            self.max_dist_c,
            self.idt_c,
            self.flag
        )
    }

    fn swap_rp(&self) -> Overlap {
        let d_left: i32;
        let d_right: i32;
        let strand1: u8;
        let bgn0: u32;
        let end0: u32;
        let bgn1: u32;
        let end1: u32;
        if self.strand1 == 0 {
            d_left = -self.d_left;
            d_right = -self.d_right;
            strand1 = 0;
            bgn0 = self.bgn0;
            end0 = self.end0;
            bgn1 = self.bgn1;
            end1 = self.end1;
        } else {
            d_left = self.d_right;
            d_right = self.d_left;
            strand1 = 1;
            bgn0 = self.len0 - self.end0;
            end0 = self.len0 - self.bgn0;
            bgn1 = self.len1 - self.end1;
            end1 = self.len1 - self.bgn1;
        }
        Overlap {
            rid0: self.rid1,
            rid1: self.rid0,
            strand1: strand1,
            len0: self.len1,
            len1: self.len0,
            d_left: d_left,
            d_right: d_right,
            bgn0: bgn1,
            end0: end1,
            bgn1: bgn0,
            end1: end0,
            dist: self.dist,
            idt: self.idt,
            dist_c: self.dist_c,
            max_dist_c: self.max_dist_c,
            idt_c: self.idt_c,
            flag: self.flag,
        }
    }

    fn reverse_strand(&self) -> Overlap {
        Overlap {
            rid0: self.rid0,
            rid1: self.rid1,
            strand1: 1 - self.strand1,
            len0: self.len1,
            len1: self.len0,
            d_left: -self.d_right,
            d_right: -self.d_left,
            bgn0: self.len0 - self.end0,
            end0: self.len0 - self.bgn0,
            bgn1: self.len1 - self.end1,
            end1: self.len1 - self.bgn1,
            dist: self.dist,
            idt: self.idt,
            dist_c: self.dist_c,
            max_dist_c: self.max_dist_c,
            idt_c: self.idt_c,
            flag: self.flag | 0x80,
        }
    }
}

type OverlapMap = FxHashMap<u32, Vec<Overlap>>;

fn build_read_ovlp_data<P>(filename: P) -> OverlapMap
where
    P: AsRef<Path>,
{
    let mut rid2ovlp = OverlapMap::default();
    let mut buffer = String::new();
    //let mut d_left: i32 = 0;

    let file = File::open(filename);
    let _err: Result<usize, io::Error> = file.unwrap().read_to_string(&mut buffer);
    for line in buffer.split("\n") {
        let mut v: Vec<&str> = Vec::<&str>::with_capacity(24); // we need pre-allocate some space for performance
        line.split(' ').for_each(|c| v.push(c));
        match v[0] {
            "O" => {
                let ovlp = Overlap::build_from(v);
                if ovlp.dist_c > 2 {
                    continue;
                }
                /*
                if ovlp.flag & 0x02 != 0x02 {
                    continue;
                }
                */
                //d_left = ovlp.d_left;
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
    let mut left_most: i32 = i32::MAX;
    let mut right_most: i32 = i32::MIN;
    let mut rlen: u32 = 0;
    for vv in v.iter() {
        if vv.dist_c > 2 {
            continue;
        }
        /*
        if vv.flag & 0x02 == 0 {
            //only consider the primary (best identity) alignments
            continue;
        }
        */
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

fn _dedup_rid2ovlap(rid2ovlp_all: FxHashMap<u32, Vec<Overlap>>) -> FxHashMap<u32, Vec<Overlap>> {
    let mut all_rids = Vec::<u32>::new();
    for r in rid2ovlp_all.keys() {
        all_rids.push(*r);
    }

    let mut best_ovlp_pair = FxHashMap::<(u32, u32), Overlap>::default();
    for r in all_rids {
        for ovlp0 in rid2ovlp_all.get(&r).unwrap() {
            let rp = (ovlp0.rid0, ovlp0.rid1);
            if best_ovlp_pair.contains_key(&rp) {
                let ovlp1 = best_ovlp_pair.get(&rp).unwrap();
                if ovlp1.dist_c < ovlp0.dist_c {
                    best_ovlp_pair.insert(rp, *ovlp0);
                }
            } else {
                best_ovlp_pair.insert(rp, *ovlp0);
            }
            // read swapped
            let rp = (ovlp0.rid1, ovlp0.rid0);
            let ovlp0_r = ovlp0.swap_rp();
            if best_ovlp_pair.contains_key(&rp) {
                let ovlp1 = best_ovlp_pair.get(&rp).unwrap();
                if ovlp1.dist_c < ovlp0_r.dist_c {
                    best_ovlp_pair.insert(rp, ovlp0_r);
                }
            } else {
                best_ovlp_pair.insert(rp, ovlp0_r);
            }
        }
    }
    let mut rid2ovlp = FxHashMap::<u32, Vec<Overlap>>::default();

    for ((rid0, _rid1), ovlp) in best_ovlp_pair {
        rid2ovlp.entry(rid0).or_insert_with(|| vec![]).push(ovlp);
    }
    rid2ovlp
}

fn get_rpair2ovlps(
    rid2ovlp_all: &FxHashMap<u32, Vec<Overlap>>,
    contained: &FxHashSet<u32>,
    chimers: &FxHashSet<u32>,
    bestn: usize,
) -> FxHashMap<ReadPair, Overlap> {
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
        for vv in v.iter() {
            if contained.contains(&vv.rid1) {
                continue;
            }
            if chimers.contains(&vv.rid1) {
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
        println!(
            "T {} {} {} {}",
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
        //let bestn = 4;
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
        if !found {
            for i in 0..bestn {
                if let Some(vv) = left_candidates.get(i) {
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
        right_candidates.sort_by(|a, b| {
            let la = a.end0 - a.bgn0;
            let lb = b.end0 - b.bgn0;
            lb.cmp(&la)
        });
        found = false;
        for i in 0..bestn {
            if let Some(vv) = right_candidates.get(i) {
                if vv.dist_c == 0 {
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
                }
            }
        }
    }
    rpair2overlap
}

fn get_g0(rpair2overlap: &FxHashMap<ReadPair, Overlap>) -> DiGraphMap<(u32, u8), u32> {
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

fn get_ulinks(
    g0: &DiGraphMap<(u32, u8), u32>,
    rpair2overlap: &FxHashMap<ReadPair, Overlap>,
) -> FxHashSet<ReadPair> {
    /*
    let mut c_nodes = FxHashSet::<(u32, u8)>::default();
    for v in g0.nodes() {
        let in_nodes = g0.neighbors_directed(v, Incoming).collect::<Vec<(u32, u8)>>();
        let out_nodes = g0.neighbors_directed(v, Outgoing).collect::<Vec<(u32, u8)>>();
        let mut count = 0;
        for v0 in in_nodes.iter() {
            for w0 in out_nodes.iter() {
                if g0.contains_edge(*v0, *w0 ) {
                    count += 1;
                    break;
                }
            }
        }
        println!("I {}:{} {}", v.0, v.1, count);
        if count == 0 {
            c_nodes.insert(v);
        }
    }
    */
    let mut u_links = FxHashSet::<ReadPair>::default();
    for n in g0.nodes() {
        /*
        if c_nodes.contains(&n) {
            continue;
        }
        */
        let mut candidates = FxHashMap::<ReadPair, Overlap>::default();
        for v in g0.neighbors_directed(n, Outgoing) {
            let rp = ReadPair {
                rid0: n.0,
                strand0: n.1,
                rid1: v.0,
                strand1: v.1,
            };
            candidates.insert(rp, *rpair2overlap.get(&rp).unwrap());
            //println!("E {} {} {} {}", n.0, n.1, v.0, v.1);
        }
        if candidates.len() == 0 {
            continue;
        }
        let mut max_ovlp_0: Option<(ReadPair, Overlap)> = None;
        let mut max_ovlp_1: Option<(ReadPair, Overlap)> = None;
        let mut max_ovlp_len_0 = 0_u32;
        let mut max_ovlp_len_1 = 0_u32;
        for (rp, vv) in candidates.iter() {
            if vv.dist_c == 0 {
                if vv.end0 - vv.bgn0 > max_ovlp_len_0 {
                    max_ovlp_len_0 = vv.end0 - vv.bgn0;
                    max_ovlp_0 = Some((*rp, *vv));
                }
            } else {
                if vv.end0 - vv.bgn0 > max_ovlp_len_1 {
                    max_ovlp_len_1 = vv.end0 - vv.bgn0;
                    max_ovlp_1 = Some((*rp, *vv));
                }
            }
        }
        if let Some((rp, _)) = max_ovlp_0 {
            u_links.insert(rp);
        } else if let Some((rp, _)) = max_ovlp_1 {
            u_links.insert(rp);
        }
    }
    u_links
}

fn get_u0(
    u_links: &FxHashSet<ReadPair>,
    _rpair2overlap: &FxHashMap<ReadPair, Overlap>,
) -> DiGraphMap<(u32, u8), ()> {
    let mut u0 = DiGraphMap::<(u32, u8), ()>::new();
    for rp in u_links.iter() {
        let rpx = rp.reverse();
        if u_links.contains(&rpx) {
            //let ovlp = rpair2overlap.get(rp).unwrap();
            //println!("C {} {}", rp.to_str(), ovlp.format());
            u0.add_edge((rp.rid0, rp.strand0), (rp.rid1, rp.strand1), ());
        }
    }
    u0
}

fn get_u0_path(u0: &DiGraphMap<(u32, u8), ()>) -> FxHashMap<u32, Vec<(u32, u8)>> {
    let mut utg0_id = 0_u32;
    let mut utg0_path = FxHashMap::<u32, Vec<(u32, u8)>>::default();
    for v in u0.nodes() {
        let in_count = u0.neighbors_directed(v, Incoming).count();
        if in_count == 0 {
            let mut path = Vec::<(u32, u8)>::new();
            let mut dfs = Dfs::new(&u0, v);
            while let Some(w) = dfs.next(&u0) {
                let out_count = u0.neighbors_directed(v, Outgoing).count();
                if out_count != 1 {
                    break;
                }
                path.push(w);
            }
            utg0_path.insert(utg0_id, path);
            utg0_id += 1;
        }
    }
    utg0_path
}

fn get_u1(
    utg0_path: &FxHashMap<u32, Vec<(u32, u8)>>,
    g0: &DiGraphMap<(u32, u8), u32>,
) -> DiGraphMap<u32, u32> {
    let mut rid2utg0 = FxHashMap::<(u32, u8), u32>::default();
    for (id, path) in utg0_path.iter() {
        for (rid, strand) in path.iter() {
            rid2utg0.insert((*rid, *strand), *id);
        }
    }

    let mut link_counter_out = FxHashMap::<u32, FxHashMap<u32, FxHashSet<u32>>>::default();
    let mut link_counter_in = FxHashMap::<u32, FxHashMap<u32, FxHashSet<u32>>>::default();
    for n in g0.nodes() {
        if !rid2utg0.contains_key(&n) {
            continue;
        }
        for v in g0.neighbors_directed(n, Outgoing) {
            if !rid2utg0.contains_key(&v) {
                continue;
            }
            let uid0 = rid2utg0.get(&n).unwrap();
            let uid1 = rid2utg0.get(&v).unwrap();
            if uid0 == uid1 {
                continue;
            }
            let c0 = link_counter_out
                .entry(*uid0)
                .or_insert(FxHashMap::<u32, FxHashSet<u32>>::default());
            let c1 = c0.entry(*uid1).or_insert(FxHashSet::<u32>::default());
            c1.insert(v.0);
            c1.insert(n.0);
        }
        for v in g0.neighbors_directed(n, Incoming) {
            if !rid2utg0.contains_key(&v) {
                continue;
            }
            let uid0 = rid2utg0.get(&n).unwrap();
            let uid1 = rid2utg0.get(&v).unwrap();
            if uid0 == uid1 {
                continue;
            }
            let c0 = link_counter_in
                .entry(*uid0)
                .or_insert(FxHashMap::<u32, FxHashSet<u32>>::default());
            let c1 = c0.entry(*uid1).or_insert(FxHashSet::<u32>::default());
            c1.insert(v.0);
            c1.insert(n.0);
        }
    }

    let mut u1 = DiGraphMap::<u32, u32>::new();
    for (uid0, c0) in link_counter_out.iter() {
        let mut max_count = 0_usize;
        for (_uid1, c1) in c0.iter() {
            if c1.len() > max_count {
                max_count = c1.len();
            }
        }
        for (uid1, c1) in c0.iter() {
            if c1.len() == max_count {
                u1.add_edge(*uid0, *uid1, max_count as u32); // reversed order for the incoming edge
            }
        }
    }

    for (uid0, c0) in link_counter_in.iter() {
        let mut max_count = 0_usize;
        for (_uid1, c1) in c0.iter() {
            if c1.len() > max_count {
                max_count = c1.len();
            }
        }
        for (uid1, c1) in c0.iter() {
            if c1.len() == max_count {
                u1.add_edge(*uid1, *uid0, max_count as u32); // reversed order for the incoming edge
            }
        }
    }
    u1
}

fn filter_u1(u1: &mut DiGraphMap<u32, u32>) -> () {
    let mut edge2remove = FxHashSet::<(u32, u32)>::default();
    for (uid0, uid1, _) in u1.all_edges() {
        if u1.contains_edge(uid1, uid0) {
            edge2remove.insert((uid0, uid1));
        }
    }
    for (uid0, uid1) in edge2remove.iter() {
        u1.remove_edge(*uid0, *uid1);
    }
    edge2remove.clear();

    for (uid0, uid1, _) in u1.all_edges() {
        let v_in_count = u1.neighbors_directed(uid0, Incoming).count();
        let v_out_count = u1.neighbors_directed(uid0, Outgoing).count();
        let w_in_count = u1.neighbors_directed(uid1, Incoming).count();
        let w_out_count = u1.neighbors_directed(uid1, Outgoing).count();
        if v_in_count >= 1 && v_out_count >= 2 && w_in_count >= 2 && w_out_count >= 1 {
            edge2remove.insert((uid0, uid1));
        }
    }
    for (uid0, uid1) in edge2remove.iter() {
        u1.remove_edge(*uid0, *uid1);
    }
}

fn bfs_extend(v: (u32, u8), g: &DiGraphMap<(u32, u8), u32>, limit: u32) -> (u32, Vec<(u32, u8)>) {
    let mut bfs = Bfs::new(g, v);
    let mut count = 0_u32;
    let mut nodes = Vec::<(u32, u8)>::with_capacity(32);
    while let Some(n) = bfs.next(g) {
        count += 1;
        if count >= limit {
            break;
        }
        nodes.push(n);
    }
    (count, nodes)
}

fn get_g1(
    u1: &DiGraphMap<u32, u32>,
    g0: &DiGraphMap<(u32, u8), u32>,
    utg0_path: &FxHashMap<u32, Vec<(u32, u8)>>,
) -> DiGraphMap<(u32, u8), u32> {
    // get G1
    let mut g1 = DiGraphMap::<(u32, u8), u32>::new();
    let mut all_edges = u1.all_edges().collect::<Vec<(u32, u32, &u32)>>();
    all_edges.sort();
    for (v, w, _) in all_edges {
        let mut out_pairs = FxHashMap::<(u32, u8), Vec<(u32, u8)>>::default();
        for vv in utg0_path.get(&v).unwrap().iter() {
            let mut max_ovlp_len = 0_u32;
            let mut max_out: Option<(u32, u8)> = None;
            for ww in utg0_path.get(&w).unwrap().iter() {
                if let Some(ovlp_len) = g0.edge_weight(*vv, *ww) {
                    if *ovlp_len > max_ovlp_len {
                        max_ovlp_len = *ovlp_len;
                        max_out = Some(*ww);
                    }
                }
            }
            if let Some(m) = max_out {
                out_pairs
                    .entry(m)
                    .or_insert(Vec::<(u32, u8)>::new())
                    .push(*vv);
            }
        }

        for (ww, vvs) in out_pairs.iter() {
            let mut max_ovlp_len = 0_u32;
            let mut max_in: Option<(u32, u8)> = None;
            for vv in vvs {
                let ovlp_len = g0.edge_weight(*vv, *ww).unwrap();
                if *ovlp_len > max_ovlp_len {
                    max_ovlp_len = *ovlp_len;
                    max_in = Some(*vv);
                }
            }
            if let Some(m) = max_in {
                let ovlp_len = g0.edge_weight(m, *ww).unwrap();
                g1.add_edge(m, *ww, *ovlp_len);
                g1.add_edge((ww.0, 1 - ww.1), (m.0, 1 - m.1), *ovlp_len);
            }
        }
    }

    for (_, path) in utg0_path {
        let mut v = path[0];
        let len = path.len();
        for w in path[1..len].iter() {
            let ovlp_len = g0.edge_weight(v, *w).unwrap();
            g1.add_edge(v, *w, *ovlp_len);
            v = *w;
        }
    }
    g1
}
fn patch_ends(
    g0: &mut DiGraphMap<(u32, u8), u32>,
    g1: &mut DiGraphMap<(u32, u8), u32>,
    _contained: &FxHashSet<u32>,
    rid2ovlp_all: &FxHashMap<u32, Vec<Overlap>>,
    rpair2overlap: &mut FxHashMap<ReadPair, Overlap>,
) -> () {
    let mut spur_nodes = Vec::<(u32, u8)>::new();
    for v in g1.nodes() {
        if g1.neighbors_directed(v, Outgoing).count() <= 1 {
            continue;
        }
        for w in g1.neighbors_directed(v, Outgoing) {
            let (count, nodes) = bfs_extend(w, &g1, 24);
            if count < 5 {
                for n in nodes {
                    spur_nodes.push(n);
                }
            }
        }
    }
    spur_nodes.sort();
    for (v0, v1) in spur_nodes {
        let v = (v0, v1);
        let mut edges2remove = FxHashSet::<((u32, u8), (u32, u8))>::default();
        for w in g1.neighbors_directed(v, Outgoing) {
            edges2remove.insert((v, w));
        }
        for w in g1.neighbors_directed(v, Incoming) {
            edges2remove.insert((w, v));
        }
        // we need to remove edges manually as there may be a bug in PetGraph that
        // the edges are not removed even the nodes are removed
        for (v, w) in edges2remove {
            g1.remove_edge(v, w);
            g1.remove_edge((w.0, 1 - w.1), (v.0, 1 - v.1));
            //println!("REMOVE {}:{} {}:{}", v.0, v.1, w.0, w.1);
            //println!("REMOVE {}:{} {}:{}", w.0, 1-w.1, v.0, 1-v.1);
        }
        let v = (v0, v1);
        g1.remove_node(v);
        //println!("REMOVE: {}:{}", v.0, v.1);
        let v = (v0, 1 - v1);
        g1.remove_node(v);
        //println!("REMOVE: {}:{}", v.0, v.1);
    }

    let mut bgn_nodes = FxHashSet::<(u32, u8)>::default();
    let mut end_nodes = FxHashSet::<(u32, u8)>::default();
    for v in g1.nodes() {
        let in_deg = g1.neighbors_directed(v, Incoming).count();
        let out_deg = g1.neighbors_directed(v, Outgoing).count();
        let mut node_set = FxHashSet::<(u32, u8)>::default();
        if in_deg == 0 && out_deg != 0 {
            //println!("DEG {}:{} {} {}", v.0, v.1, in_deg, out_deg);
            bgn_nodes.insert(v);
            end_nodes.insert((v.0, 1 - v.1));
            // println!("YB {}:{}", v.0, v.1);
            // println!("YE {}:{}", v.0, 1 - v.1);
            let (_count, nodes) = bfs_extend(v, &g1, 8);
            for n in nodes {
                node_set.insert(n);
            }
            let mut node_vec = node_set.into_iter().collect::<Vec<(u32, u8)>>();
            node_vec.sort();
            for n in node_vec {
                let mut flag = false;
                for v in g0.neighbors_directed(n, Incoming) {
                    if !g1.contains_edge(v, n) {
                        flag = true;
                        break;
                    }
                }
                if flag {
                    bgn_nodes.insert(n);
                    end_nodes.insert((n.0, 1 - n.1));
                    // println!("YB {}:{}", n.0, n.1);
                    // println!("YE {}:{}", n.0, 1 - n.1);
                }
            }
        }
    }

    // handle some false containment cases
    //let mut g_tmp = DiGraphMap::<(u32, u8), u32>::new();
    //let mut m_read_ids = HashSet::<u32>::new();
    let mut end_nodes_vec = end_nodes.into_iter().collect::<Vec<(u32, u8)>>();
    end_nodes_vec.sort();
    for r in end_nodes_vec.iter() {
        if !rid2ovlp_all.contains_key(&r.0) {
            continue;
        }

        // if g0.neighbors_directed(*r, Outgoing).count() > 0 {
        //    continue;
        // }

        let v = rid2ovlp_all.get(&r.0).unwrap();
        for vv in v.iter() {
            let rid0 = vv.rid0;
            let strand0 = 0;
            let rid1 = vv.rid1;
            let strand1 = vv.strand1;
            let ovlp_len = vv.end0 - vv.bgn0;
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

    fn find_path(
        g: &DiGraphMap<(u32, u8), u32>,
        v: &(u32, u8),
        w: &(u32, u8),
        limit: u32,
    ) -> Option<Vec<(u32, u8)>> {
        let mut bfs = Bfs::new(&g, (v.0, v.1));
        let mut max_depth = FxHashMap::<(u32, u8), u32>::default();
        let mut max_pre = FxHashMap::<(u32, u8), (u32, u8)>::default();
        let mut count = 0_u32;
        let mut found = false;
        while let Some(n) = bfs.next(&g) {
            let mut m = 0_u32;
            let mut p = (0_u32, 255_u8);
            for ww in g.neighbors_directed(n, Incoming) {
                if let Some(m0) = max_depth.get(&ww) {
                    if *m0 > m {
                        m = *m0;
                        p = ww;
                    }
                }
            }
            max_depth.insert(n, m + 1);
            max_pre.insert(n, p);
            if n == *w {
                found = true;
                break;
            }
            count += 1;
            if count > limit {
                break;
            }
        }
        if found {
            let mut path = Vec::<(u32, u8)>::new();
            let mut vv = w;
            while vv != v && vv.1 != 255 {
                path.push(*vv);
                vv = max_pre.get(vv).unwrap();
            }
            path.push(*vv);
            path.reverse();
            Some(path)
        } else {
            None
        }
    }

    for v in end_nodes_vec {
        let (mut _count, mut nodes) = bfs_extend(v, g0, 48);
        nodes.reverse();
        for w in nodes {
            // println!("try: {}:{} {}:{}", v.0, v.1, w.0, w.1);
            if bgn_nodes.contains(&w) {
                if let Some(_path) = find_path(&g1, &v, &w, 48) {
                    // println!("try: {}:{} {}:{}, path exists in g1", v.0, v.1, w.0, w.1);
                    // for vv in path {
                    //     println!("try: {}:{} {}:{}, p:{}:{}", v.0, v.1, w.0, w.1, vv.0, vv.1);
                    // }
                    continue;
                }
                if let Some(path) = find_path(g0, &v, &w, 48) {
                    let mut vv = path[0];
                    let vlen = path.len();
                    for ww in path[1..vlen].iter() {
                        if !g1.contains_edge(v, *ww) {
                            let ovlp_len = g0.edge_weight(vv, *ww).unwrap();
                            g1.add_edge(vv, *ww, *ovlp_len);
                            g1.add_edge((ww.0, 1 - ww.1), (vv.0, 1 - vv.1), *ovlp_len);
                            // println!("ADD0 {}:{} {}:{}", vv.0, vv.1, ww.0, ww.1);
                            // println!("ADD1 {}:{} {}:{}", ww.0, 1 - ww.1, vv.0, 1 - vv.1);
                        }
                        vv = *ww;
                    }
                    break;
                };
            }
        }
    }
}

fn find_branching_nodes(
    g: &mut DiGraphMap<(u32, u8), u32>,
    candidates: Vec<(u32, u8)>,
    max_path_length: u32,
    min_path_length: u32,
) -> Vec<(u32, u8)> {
    let mut branching_nodes = FxHashSet::<(u32, u8)>::default();
    for v in candidates {
        let out_count = g.neighbors_directed(v, Outgoing).count();
        if out_count >= 2 {
            let mut overlapped_path = false;
            let mut ext_branch_count = 0_u32;
            let mut node_count = FxHashMap::<(u32, u8), u32>::default();
            for w in g.neighbors_directed(v, Outgoing) {
                let (count, p) = bfs_extend(w, g, max_path_length);
                if count > min_path_length {
                    ext_branch_count += 1;
                    for vv in p {
                        let c = node_count.entry(vv).or_insert(0);
                        *c += 1;
                        if *c > 1 {
                            overlapped_path = true;
                            break;
                        }
                    }
                }
                if overlapped_path {
                    break;
                }
            }
            if ext_branch_count >= 2 && !overlapped_path {
                branching_nodes.insert(v);
            }
        }
    }
    branching_nodes.iter().map(|x| *x).collect()
}

fn bfs_shortest(
    v: (u32, u8),
    g: &DiGraphMap<(u32, u8), u32>,
    limit: u32,
) -> (FxHashMap<(u32, u8), u32>, FxHashMap<(u32, u8), (u32, u8)>) {
    let mut bfs = Bfs::new(g, v);
    //let mut count = 0_u32;
    let mut lengths = FxHashMap::<(u32, u8), u32>::default();
    let mut bt_node = FxHashMap::<(u32, u8), (u32, u8)>::default();
    lengths.insert(v, 0);
    while let Some(n) = bfs.next(g) {
        if n == v {
            continue;
        }
        let mut min_length = std::u32::MAX;
        let mut min_p_node: Option<(u32, u8)> = None;
        for w in g.neighbors_directed(n, Incoming) {
            if let Some(l) = lengths.get(&w) {
                if *l < min_length {
                    min_p_node = Some(w);
                    min_length = *l;
                }
            }
        }
        if min_length > limit {
            break;
        }
        if let Some(p_node) = min_p_node {
            bt_node.insert(n, p_node);
            lengths.insert(n, min_length + 1);
        }
    }
    (lengths, bt_node)
}

fn remove_repeat_branch(g: &mut DiGraphMap<(u32, u8), u32>) -> () {
    let candidates = g.nodes().collect();
    let mut branching_nodes = find_branching_nodes(g, candidates, 64, 56);
    branching_nodes.sort();
    for v in branching_nodes {
        let mut best_ovlp_len = 0_u32;
        let mut best_node: Option<(u32, u8)> = None;
        let mut out_nodes = g
            .neighbors_directed(v, Outgoing)
            .collect::<Vec<(u32, u8)>>();
        out_nodes.sort();
        for w in out_nodes {
            let ovlp_len = *g.edge_weight(v, w).unwrap();
            if ovlp_len > best_ovlp_len {
                best_ovlp_len = ovlp_len;
                best_node = Some(w);
            }
        }
        let mut edges2remove = FxHashSet::<((u32, u8), (u32, u8))>::default();
        if let Some(ww) = best_node {
            for w in g.neighbors_directed(v, Outgoing) {
                if w != ww {
                    edges2remove.insert((v, w));
                    edges2remove.insert(((w.0, 1 - w.1), (v.0, 1 - v.1)));
                }
            }
        }
        for (v, w) in edges2remove {
            g.remove_edge(v, w);
        }
    }
}

fn remove_repeat_bridge(
    g: &mut DiGraphMap<(u32, u8), u32>,
    candidates: Vec<(u32, u8)>,
) -> Vec<(u32, u8)> {
    let mut branching_nodes = find_branching_nodes(g, candidates, 32, 28);
    let mut branching_nodes_set = FxHashSet::<u32>::default();
    for v in branching_nodes.iter() {
        branching_nodes_set.insert(v.0);
    }
    branching_nodes.sort();
    let mut new_candidates = FxHashSet::<(u32, u8)>::default();
    for v in branching_nodes {
        // println!("V {}:{}", v.0, v.1);
        let (lengths, bt_node) = bfs_shortest(v, &g, 24);
        let mut len_vec = Vec::<(u32, u32, u8)>::new();
        for (w, len) in lengths {
            if w == v {
                continue;
            }
            // let out_count = g.neighbors_directed(w, Outgoing).count();
            // let in_count = g.neighbors_directed(w, Incoming).count();
            // println!("W {}:{} {} {} {}", w.0, w.1, len, in_count, out_count);
            if branching_nodes_set.contains(&w.0) {
                len_vec.push((len, w.0, w.1));
                // println!("W2 {}:{} {} {} {}", w.0, w.1, len, in_count, out_count);
            }
        }
        if len_vec.len() == 0 {
            continue;
        }
        len_vec.sort();
        // println!("WL {}", len_vec.len());
        let mut path = Vec::<(u32, u8)>::new();
        let mut vv = (len_vec[0].1, len_vec[0].2);
        path.push(vv);
        loop {
            if let Some(ww) = bt_node.get(&vv) {
                path.push(*ww);
                vv = *ww;
            } else {
                break;
            }
        }
        let plen = path.len();
        // println!("P {}", plen);
        if plen > 8 {
            continue;
        }
        path.reverse();
        for vv in path.iter() {
            for w in g.neighbors(*vv) {
                new_candidates.insert(w);
            }
            for w in g.neighbors((vv.0, 1 - vv.1)) {
                new_candidates.insert(w);
            }
        }

        vv = path[0];

        for ww in path[1..plen].iter() {
            // println!("R {}:{} {}:{}", vv.0, vv.1, ww.0, ww.1);
            g.remove_edge(vv, *ww);
            g.remove_edge((ww.0, 1 - ww.1), (vv.0, 1 - vv.1));
            vv = *ww;
        }
    }
    new_candidates.iter().map(|x| *x).collect()
}

fn get_path_from_seed(
    g1: &DiGraphMap<(u32, u8), u32>,
    seed: &(u32, u8, u32),
    used_reads: &mut FxHashSet<u32>,
) -> Vec<((u32, u8), (u32, u8))> {
    let v = *seed;
    let mut bfs = Bfs::new(&g1, (v.0, v.1));
    let mut max_depth = FxHashMap::<(u32, u8), u32>::default();
    let mut max_pre = FxHashMap::<(u32, u8), (u32, u8)>::default();
    let mut end_nodes = Vec::<(u32, u32, u8)>::new();
    max_depth.insert((v.0, v.1), 0);

    let mut mm = 0_u32;
    let mut nn: Option<(u32, u8)> = None;
    loop {
        // we will eventully need a better graph travesal algorithm
        if let Some(n) = bfs.next(&g1) {
            let mut m = 0_u32;
            let mut p = (0_u32, 255_u8);
            for w in g1.neighbors_directed(n, Incoming) {
                if let Some(m0) = max_depth.get(&w) {
                    if *m0 > m {
                        m = *m0;
                        p = w;
                    }
                }
            }
            nn = Some((n.0, n.1));
            mm = m;
            max_depth.insert(n, m + 1);
            max_pre.insert(n, p);

            let out_count = g1.neighbors_directed((n.0, n.1), Outgoing).count();
            if out_count == 0 {
                end_nodes.push((m + 1, n.0, n.1));
            } else if used_reads.contains(&n.0) {
                end_nodes.push((m + 1, n.0, n.1));
                break;
            }
        } else {
            // traversal end
            if let Some(n) = nn {
                end_nodes.push((mm + 1, n.0, n.1));
            }
            break;
        }
    }

    let mut edge_list = Vec::<((u32, u8), (u32, u8))>::new();
    if end_nodes.len() == 0 {
        return edge_list;
    }
    end_nodes.sort_by(|a, b| (b.0).partial_cmp(&a.0).unwrap());
    let mut w = (end_nodes[0].1, end_nodes[0].2);
    if end_nodes[0].0 < 10 {
        return edge_list;
    }

    while w.1 != 255 {
        let v = max_pre.get(&w).unwrap();
        edge_list.push(((v.0, v.1), (w.0, w.1)));
        w = *v;
    }
    edge_list.reverse();
    edge_list
}

fn write_paths_from_seeds(
    g1: &DiGraphMap<(u32, u8), u32>,
    seeds: &Vec<(u32, u8, u32)>,
    mut used_reads: &mut FxHashSet<u32>,
    layout_file: &mut BufWriter<File>,
    rpair2overlap: &FxHashMap<ReadPair,Overlap>,
    utg_id: &mut u32,
) -> () {
    for v in seeds {
        if used_reads.contains(&v.0) {
            continue;
        }
        let edge_list = get_path_from_seed(&g1, &v, &mut used_reads);
        if edge_list.len() == 0 {
            continue;
        }
        let mut edge_list2 = Vec::<((u32, u8), (u32, u8))>::with_capacity(128); 
        let mut used_reads2 = FxHashSet::<u32>::default();
        for (v0, v1) in edge_list {
            if used_reads2.contains(&v0.0) {
                break;
            }
            edge_list2.push((v0, v1));
            used_reads2.insert(v0.0);
        }
        let (_, w) = edge_list2[edge_list2.len() - 1];
        let _res = writeln!(
            layout_file,
            "U {} {}:{} {} {}:{} {}",
            utg_id, v.0, v.1, v.2, w.0, w.1, edge_list2.len()
        );
        let mut used_reads_in_path = FxHashSet::<u32>::default();
        for (v, w) in edge_list2 {
            if v.1 == 255 {
                continue;
            }
            if used_reads_in_path.contains(&v.0) {
                break;
            }
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
                "E {} {} {} {} {} {} {}",
                utg_id, v.0, v.1, w.0, w.1, ovlp.bgn0, ovlp.bgn1
            );
            used_reads_in_path.insert(v.0);
            used_reads.insert(w.0);
        }
        *utg_id += 1;
    }
}

fn generate_layout(
    g1: &DiGraphMap<(u32, u8), u32>,
    rpair2overlap: &FxHashMap<ReadPair, Overlap>,
    output_prefix: &str,
) -> () {
    let mut seeds = Vec::<(u32, u8, u32)>::new(); // a map from each component to the seed nodes
    for v in g1.nodes() {
        let in_degree = g1.neighbors_directed(v, Incoming).count();
        let out_degree = g1.neighbors_directed(v, Outgoing).count();
        if in_degree == 0 && out_degree != 0 {
            let (count, _nodes) = bfs_extend(v, &g1, 1000000);
            seeds.push((v.0, v.1, count));
        }
    }

    let mut used_reads = FxHashSet::<u32>::default();
    seeds.sort_by(|a, b| (b.2).partial_cmp(&a.2).unwrap());
    let layout_filename = format!("{}_layout.dat", output_prefix);
    let mut layout_file = BufWriter::new(File::create(layout_filename).unwrap());
    let mut utg_id = 0_u32;
    write_paths_from_seeds(&g1, &seeds, &mut used_reads, &mut layout_file, &rpair2overlap, &mut utg_id);

    // process the rest un-layouted nodes, likely in circles

    let mut g1_res = DiGraphMap::<(u32, u8), u32>::new();
    for e in g1.all_edges() {
        if used_reads.contains(&e.0.0) || used_reads.contains(&e.1.0) {
            continue;
        } 
        g1_res.add_edge(e.0, e.1, *e.2);
    }

    let mut seeds = Vec::<(u32, u8, u32)>::new(); // a map from each component to the seed nodes
    let mut reachable = FxHashSet::<(u32, u8)>::default();
    for v in g1_res.nodes() {
        if used_reads.contains(&v.0) {
            continue;
        }
        if reachable.contains(&v) {
            continue;
        }
        let (count, _nodes) = bfs_extend(v, &g1_res, 1000000);
        for w in _nodes {
            reachable.insert(w);
        }
        seeds.push((v.0, v.1, count));
    }
    seeds.sort_by(|a, b| (b.2).partial_cmp(&a.2).unwrap());
    write_paths_from_seeds(&g1_res, &seeds, &mut used_reads, &mut layout_file, &rpair2overlap, &mut utg_id); 
}

pub fn ovlp2layout_v1(prefix: &String, out_prefix: &String, bestn: usize) -> () {
    let prefix = prefix.clone();
    let infile_pattern = [prefix, "*".to_string()].concat();

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
    //let mut rid2ovlp = OverlapMap::new();
    let mut rid2ovlp_all = OverlapMap::default();
    for child in children {
        let rid2ovlp_p = child.join().expect("oops! the child thread panicked");
        rid2ovlp_all.extend(rid2ovlp_p);
    }
    //let rid2ovlp_all = dedup_rid2ovlap(rid2ovlp_all);

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

    /*
    println!(
        "#LEN {} {} {} {}",
        rid2ovlp_all.len(),
        contained.len(),
        chimers.len(),
        dead_ends.len()
    );
    */

    /*
    let dead_end_file = File::create("./dead_ends").unwrap();

    for rid in dead_ends.iter() {
        let _err = writeln!(&dead_end_file, "{}", rid);
    }
    */

    // graph processing
    let mut rpair2overlap = get_rpair2ovlps(&rid2ovlp_all, &contained, &chimers, bestn);
    let mut g0 = get_g0(&rpair2overlap);
    let u_links = get_ulinks(&g0, &rpair2overlap);
    let u0 = get_u0(&u_links, &rpair2overlap);
    let utg0_path = get_u0_path(&u0);
    let mut u1 = get_u1(&utg0_path, &g0);
    filter_u1(&mut u1);

    /*
    for (id, path) in utg0_path.iter() {
        print!("P {}", id);
        for (rid, strand) in path.iter() {
            print!(" {} {}", rid, strand);
        }
        println!();
    }

    for (uid1, uid0, _) in u1.all_edges() {
        println!("L {} {}", uid1, uid0);
    }
    */

    // expanding utgs
    let mut g1 = get_g1(&u1, &g0, &utg0_path);

    // println!("# {} {}", g0.nodes().count(), g1.nodes().count());

    let candidates: Vec<(u32, u8)> = g1.nodes().collect();
    remove_repeat_bridge(&mut g1, candidates);
    remove_repeat_branch(&mut g1);
    patch_ends(
        &mut g0,
        &mut g1,
        &contained,
        &rid2ovlp_all,
        &mut rpair2overlap,
    );
    let candidates: Vec<(u32, u8)> = g1.nodes().collect();
    let candidates = remove_repeat_bridge(&mut g1, candidates);
    let candidates = remove_repeat_bridge(&mut g1, candidates);
    remove_repeat_bridge(&mut g1, candidates);
    remove_repeat_branch(&mut g1);

    let g0_filename = format!("{}_g0.dat", out_prefix);
    let mut graph0_file = BufWriter::new(File::create(g0_filename).unwrap());
    for (v, w, len) in g0.all_edges() {
        let _res = writeln!(graph0_file, "G {} {} {} {} {}", v.0, v.1, w.0, w.1, len);
    }

    let g1_filename = format!("{}_g1.dat", out_prefix);
    let mut graph1_file = BufWriter::new(File::create(g1_filename).unwrap());
    // let mut pair_file = BufWriter::new(File::create("pair1.dat").unwrap());
    for (v, w, len) in g1.all_edges() {
        let _res = writeln!(graph1_file, "G {} {} {} {} {}", v.0, v.1, w.0, w.1, len);
        /*
        let rp = ReadPair {
            rid0: v.0,
            strand0: v.1,
            rid1: w.0,
            strand1: w.1,
        };

        let ovlp = *rpair2overlap.get(&rp).unwrap();
        let _res = writeln!(
            pair_file,
            "G {} {} {} {} {} {}",
            v.0,
            v.1,
            w.0,
            w.1,
            len,
            ovlp.format()
        );
        */
    }
    generate_layout(&g1, &rpair2overlap, &out_prefix);
}
