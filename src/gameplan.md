First, we need to generate all the nets.
We'll take the kth net, and find the flatness of each point.
We'll make a vector of the flat points (alpha < 1/16).
As we do this, we should make a map<pair<int,point<2> (or Point)>,cylinder>.
Essentially, we're making a structure that maps a particular point in a net to its cylinder.
The cylinder stores the flatness in it, as well as two points.

Since we calculated the cylinders, we've got all the lines.

calc nets


for (net : nets) {

}
