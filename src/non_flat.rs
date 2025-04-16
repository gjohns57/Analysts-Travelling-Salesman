use std::collections::{HashSet, VecDeque};
use nalgebra::Vector2;
use crate::graph::{AdjacencyMatrix, Graph};

const C0: f32 = 0.2;

/*fn get_points_in_ball(points: &Vec<Vector2<f32>>, net: &Vec<usize>, center: usize, radius: f32) -> HashSet<usize> {
    net.iter().map(|index_ref| *index_ref).filter(|index| (points[*index] - points[center]).norm_squared() < radius * radius).collect()
}*/

fn get_points_in_ball(_points: &Vec<Vector2<f32>>, _net: &Vec<usize>, _center: usize, _radius: f32) -> HashSet<usize> {
    let mut temp: HashSet<usize> = HashSet::new();

    temp.insert(0);
    temp.insert(1);
    temp.insert(2);
    temp.insert(3);
    temp.insert(4);
    temp.insert(5);
    temp.insert(6);
    temp.insert(7);

    temp
}

fn bfs(
    start: usize, 
    visited: &mut Vec<bool>, 
    reps: &mut Vec<usize>, 
    graph: &mut AdjacencyMatrix, 
    used: &mut AdjacencyMatrix,
    origin: usize,
) {
    let mut q: VecDeque<usize> = VecDeque::new();
    q.push_back(start);

    while !q.is_empty() {
        let node = *q.front().unwrap();
        
        //Pops out the rep if the origin is inside the component
        //So we dont connect origin to its own component
        if node == origin {
            reps.pop();
        }

        //I think flat pairs should be just a line, but
        //kept this just in case
        if visited[node] {
            continue;
        }

        visited[node] = true;
        q.pop_front();

        for i in  0..graph.vertex_ct() {
            if graph.adjacent(node, i) && !visited[i] {
                //Checks if the edge is in used
                if used.adjacent(node, i) {
                    //Removes it from the graph if it is
                    graph.remove_edge(node, i);
                }
                else {
                    //Adds the edge to the used and pushes back the vertex
                    used.add_edge(node, i);
                    q.push_back(i);
                }
            }
        }
    }
}

fn find_reps(edges: &mut AdjacencyMatrix, used: &mut AdjacencyMatrix, origin: usize) -> Vec<usize> {
    let mut reps: Vec<usize> = Vec::new();
    let mut visited = vec![false; edges.vertex_ct()];

    //Goes through each vertex and finds a vertex from each component
    //Will remove edges already used
    for i in 0..edges.vertex_ct() {
        if visited[i] {
            continue;
        }
        reps.push(i);
        bfs(i, &mut visited, &mut reps, edges, used, origin);
    }

    reps
}

fn connect(origin: usize, reps: &Vec<usize>, edges: &mut AdjacencyMatrix, used: &mut AdjacencyMatrix) {
    //Connects the origin to each component rep
    for i in 0..reps.len() {
        edges.add_edge(origin, reps[i]);
        used.add_edge(origin, reps[i]);
    }
}

fn has_edge(used: &AdjacencyMatrix, v: usize) -> bool {
    for i in 0..used.vertex_ct() {
        if used.adjacent(i, v) {
            return true
        }
    }

    false
}

fn contains_edge(graph: &AdjacencyMatrix, verts: &HashSet<usize>) -> HashSet<usize> {
    let mut contains: HashSet<usize> = HashSet::new();

    for v in verts {
        if has_edge(graph, *v) {
            contains.insert(*v);
        }
    }

    contains
}

//This function should add edges between flat pairs
fn add_flat_pairs(
    graph: &mut AdjacencyMatrix, 
    _node: usize, 
    _net1: &Vec<usize>, 
    _net2: &Vec<usize>,
    verts: &mut HashSet<usize>,
) {
    //For testing we will say flat pairs were 1/2 and 3/4
    graph.add_edge(1, 2);
    graph.add_edge(3,4);
    graph.add_edge(4, 7);
    verts.insert(1);
    verts.insert(2);
    verts.insert(3);
    verts.insert(4);
    verts.insert(7);
}

pub fn non_flat(
    points: &Vec<Vector2<f32>>,
    graph: &mut AdjacencyMatrix, //next_graph
    net_k: &Vec<usize>,
    net_k_plus_1: &Vec<usize>,
    not_flat_k: &Vec<usize>,
    _n_k: i32,
    n_k_plus_1: i32,
    r0: f32,
) {

    for u in not_flat_k {
        //Index of all points in net_k_plus_1 that are inside the ball of the current u
        let vp_k_plus_1 = get_points_in_ball(points, net_k_plus_1, *u, C0 * 2.0f32.powi(-n_k_plus_1) * r0);

        //Contains edge should return a hashset of the indices that have an edge in next_graph
        //So we take the difference - the contains_edge function might need to be rewritten, I wasnt sure
        //exactly what to do with the AdjacencyMatrix
        let v_k_plus_1: HashSet<usize> = vp_k_plus_1.difference(&contains_edge(graph, &vp_k_plus_1)).cloned().collect();

        for t in &v_k_plus_1 {
            print!("{}", *t);
        }
        println!();

        //Goes to next point if its empty
        if v_k_plus_1.is_empty() {
            continue;
        }

        let mut g_k_plus_1 = AdjacencyMatrix::new();
        g_k_plus_1.resize(net_k_plus_1.len());
        let mut in_flat_pair: HashSet<usize> = HashSet::new();

        for v in v_k_plus_1.iter() {
            add_flat_pairs(&mut g_k_plus_1, *v, net_k, net_k_plus_1, &mut in_flat_pair);
        }

        //Only issue is new_graph is size net_k_plus_1 and not vs_k_plus_1
        //But I had to resize it to add flat pairs to it, but I need flat pairs to make vs
        //Somewhat circular, maybe I'm being dumb
        //let vs_k_plus_1: HashSet<usize> = vp_k_plus_1.union(&in_flat_pair).cloned().collect();

        let mut reps = find_reps(&mut g_k_plus_1, graph, *u);
        connect(*u, &mut reps, &mut g_k_plus_1, graph);
        g_k_plus_1.print_matrix();
    }
}