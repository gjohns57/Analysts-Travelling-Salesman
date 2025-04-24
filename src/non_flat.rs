use std::collections::{HashSet, VecDeque};
use nalgebra::Vector2;
use crate::graph::{AdjacencyMatrix, Graph};

const C0: f32 = 0.2;

#[derive(Default)]
struct Cylinder {
    x0: Vector2<f32>,
    x1: Vector2<f32>,
    width: f32,
}

//Can remove this, I was having issues connecting
fn find_thinnest_cylinder(points: &Vec<Vector2<f32>>, net: &Vec<usize>, center: usize, radius: f32) -> Cylinder {
    let side_offsets = vec![radius * Vector2::<f32>::x(), radius * Vector2::<f32>::y(), -radius * Vector2::<f32>::x(), -radius * Vector2::<f32>::y()];
    let l = (40. * C0).ceil() as i32;
    let mut cylinder_width = f32::MAX;
    let points_in_ball = get_points_in_ball(points, net, center, radius);
    let mut cylinder = Cylinder::default();

    // Atrocious nesting
    for side1 in 0..side_offsets.len() {
        for side2 in 0..side1 {
            for i in (-2 * l + 1)..(2 * l) {
                let x0 = side_offsets[side1] + Vector2::<f32>::new(side_offsets[side1].y, -side_offsets[side1].x) * radius * ((i as f32) / ((2 * l) as f32));
                for j in (-2 * l + 1)..(2 * l) {
                    let x1 = side_offsets[side2] + Vector2::<f32>::new(side_offsets[side2].y, -side_offsets[side2].x) * radius * ((j as f32) / ((2 * l) as f32));

                    let distance = max_distance_from_line(points, x0, x1, &points_in_ball);

                    if distance > cylinder_width {
                        cylinder_width = distance;
                        cylinder.width = cylinder_width;
                        cylinder.x0 = x0;
                        cylinder.x1 = x1;
                    }
                }
            }
        }
    }

    cylinder
}

//Can remove this, I was just having issues connecting
fn get_points_in_ball(points: &Vec<Vector2<f32>>, net: &Vec<usize>, center: usize, radius: f32) -> HashSet<usize> {
    net.iter().map(|index_ref| *index_ref).filter(|index| (points[*index] - points[center]).norm_squared() < radius * radius).collect()
}

fn bfs(
    start: usize, 
    visited: &mut HashSet<usize>, 
    reps: &mut Vec<usize>, 
    graph: &mut AdjacencyMatrix, 
    used: &mut AdjacencyMatrix,
    set: &HashSet<usize>,
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
        if visited.contains(&node) {
            continue;
        }

        visited.insert(node);
        q.pop_front();

        for i in set {
            if graph.adjacent(node, *i) && !visited.contains(i) {
                //Checks if the edge is in used
                if used.adjacent(node, *i) {
                    //Removes it from the graph if it is
                    graph.remove_edge(node, *i);
                }
                else {
                    //Adds the edge to the used and pushes back the vertex
                    used.add_edge(node, *i);
                    q.push_back(*i);
                }
            }
        }
    }
}

fn find_reps(edges: &mut AdjacencyMatrix, used: &mut AdjacencyMatrix, set: &HashSet<usize>, origin: usize) -> Vec<usize> {
    let mut reps: Vec<usize> = Vec::new();
    let mut visited: HashSet<usize> = HashSet::new();

    //Goes through each vertex and finds a vertex from each component
    //Will remove edges already used
    for i in set {
        if visited.contains(i) {
            continue;
        }
        reps.push(i);
        bfs(*i, &mut visited, &mut reps, edges, used, &set, origin);
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
    node: usize, 
    points: &Vec<Vector2<f32>>,
    net1: &Vec<usize>, 
    net2: &Vec<usize>,
    verts: &mut HashSet<usize>,
    epsilon: f32,
    k: i32,
) {
    //Also, feels like I should only be iterating
    //in here or in the non_flat function not both
    //Might just want to make center the node and iterate in non_flat
    for p in *net1 {
        let xp_in_ball: HashSet<usize> = get_points_in_ball(points, net2, p, epsilon);
        
        //This feels bad but thinnest cylinder takes a vec
        let xp_in_ball_vec: Vec<usize> = xp_in_ball.into_iter().collect();

        let c: Cylinder = find_thinnest_cylinder(points, &xp_in_ball_vec, p, epsilon);
        let alpha: f64 = (1 << k) * c.width / epsilon;

        if alpha < 1.0 / 16 {
            let mut min_cylinder_aligned_component: f64 = 1.0e300;
            
            //Set this to p so that if no flat pair then the first part of the && will fail
            let mut next_point: usize = p;
            for q in xp_in_ball {
                if (points[p] - points[q]).norm() < epsilon || (points[q] - points[p]).norm() >= C0 * powi(2, -k - 1) * epsilon {
                    continue;
                }

                let component: f64 = ((points[q] - points[p]) * (c.x1 - c.x0) / (c.x1 - c.x0).norm()).abs();
                if component < min_cylinder_aligned_component {
                    min_cylinder_aligned_component = component;
                    next_point = q;
                }
            }

            if (points[next_point] - points[p]).norm() > epsilon && (points[next_point] - points[p]).norm() < 2 * epsilon {
                graph.add_edge(p, next_point);
                verts.insert(p);
                verts.insert(next_point);
            }
        }
    }
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
        //So we take the difference
        let v_k_plus_1: HashSet<usize> = vp_k_plus_1.difference(&contains_edge(graph, &vp_k_plus_1)).cloned().collect();

        //Goes to next point if its empty
        if v_k_plus_1.is_empty() {
            continue;
        }

        let mut g_k_plus_1 = AdjacencyMatrix::new();
        g_k_plus_1.resize(net_k_plus_1.len());
        let mut in_flat_pair: HashSet<usize> = HashSet::new();

        /* I was not sure what exactly to pass into epsilon and k */
        for v in v_k_plus_1.iter() {
            add_flat_pairs(&mut g_k_plus_1, *v, &points, net_k, net_k_plus_1, &mut in_flat_pair, r0, n_k_plus_1);
        }

        let vs_k_plus_1: HashSet<usize> = vp_k_plus_1.union(&in_flat_pair).cloned().collect();

        let mut reps = find_reps(&mut g_k_plus_1, graph, &vs_k_plus_1, *u);
        connect(*u, &mut reps, &mut g_k_plus_1, graph);
        g_k_plus_1.print_matrix();
    }
}