import {addBeach, removeBeach} from "./Beach";
import {sortCellHalfedges, cellHalfedgeStart, clipCells} from "./Cell";
import {firstCircle} from "./Circle";
import {clipEdges} from "./Edge";
import RedBlackTree from "./RedBlackTree";

export var epsilon = 1e-6;
export var epsilon2 = 1e-12;
export var beaches;
export var cells;
export var circles;
export var edges;

function triangleArea(a, b, c) {
  return (a[0] - c[0]) * (b[1] - a[1]) - (a[0] - b[0]) * (c[1] - a[1]);
}

function lexicographic(a, b) {
  return b[1] - a[1]
      || b[0] - a[0];
}

export default function Diagram(sites, extent) {
  var site = sites.sort(lexicographic).pop(),
      x,
      y,
      circle;

  edges = [];
  cells = new Array(sites.length);
  beaches = new RedBlackTree;
  circles = new RedBlackTree;

  while (true) {
    circle = firstCircle;
    if (site && (!circle || site[1] < circle.y || (site[1] === circle.y && site[0] < circle.x))) {
      if (site[0] !== x || site[1] !== y) {
        addBeach(site);
        x = site[0], y = site[1];
      }
      site = sites.pop();
    } else if (circle) {
      removeBeach(circle.arc);
    } else {
      break;
    }
  }

  sortCellHalfedges();

  if (extent) {
    var x0 = +extent[0][0],
        y0 = +extent[0][1],
        x1 = +extent[1][0],
        y1 = +extent[1][1];
    clipEdges(x0, y0, x1, y1);
    clipCells(x0, y0, x1, y1);
  }

  this.edges = edges;
  this.cells = cells;

  this.delaunayTriangles =
  beaches =
  circles =
  edges =
  cells = null;
}

function supertriangle(vertices) {
  var xmin = Number.POSITIVE_INFINITY,
      ymin = Number.POSITIVE_INFINITY,
      xmax = Number.NEGATIVE_INFINITY,
      ymax = Number.NEGATIVE_INFINITY,
      i, dx, dy, dmax, xmid, ymid;

  for (i = vertices.length; i--;) {
    if (vertices[i][0] < xmin) xmin = vertices[i][0];
    if (vertices[i][0] > xmax) xmax = vertices[i][0];
    if (vertices[i][1] < ymin) ymin = vertices[i][1];
    if (vertices[i][1] > ymax) ymax = vertices[i][1];
  }

  dx = xmax - xmin;
  dy = ymax - ymin;
  dmax = Math.max(dx, dy);
  xmid = xmin + dx * 0.5;
  ymid = ymin + dy * 0.5;

  return [
    [xmid - 20 * dmax, ymid - dmax],
    [xmid, ymid + 20 * dmax],
    [xmid + 20 * dmax, ymid - dmax]
  ];
}

function circumcircle(vertices, i, j, k) {
  var x1 = vertices[i][0],
      y1 = vertices[i][1],
      x2 = vertices[j][0],
      y2 = vertices[j][1],
      x3 = vertices[k][0],
      y3 = vertices[k][1],
      fabsy1y2 = Math.abs(y1 - y2),
      fabsy2y3 = Math.abs(y2 - y3),
      xc, yc, m1, m2, mx1, mx2, my1, my2, dx, dy;

  /* Check for coincident points */
  if (fabsy1y2 < epsilon && fabsy2y3 < epsilon)
    throw new Error("Eek! Coincident points!");

  if (fabsy1y2 < epsilon) {
    m2 = -((x3 - x2) / (y3 - y2));
    mx2 = (x2 + x3) / 2.0;
    my2 = (y2 + y3) / 2.0;
    xc = (x2 + x1) / 2.0;
    yc = m2 * (xc - mx2) + my2;
  }

  else if (fabsy2y3 < epsilon) {
    m1 = -((x2 - x1) / (y2 - y1));
    mx1 = (x1 + x2) / 2.0;
    my1 = (y1 + y2) / 2.0;
    xc = (x3 + x2) / 2.0;
    yc = m1 * (xc - mx1) + my1;
  }

  else {
    m1 = -((x2 - x1) / (y2 - y1));
    m2 = -((x3 - x2) / (y3 - y2));
    mx1 = (x1 + x2) / 2.0;
    mx2 = (x2 + x3) / 2.0;
    my1 = (y1 + y2) / 2.0;
    my2 = (y2 + y3) / 2.0;
    xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
    yc = (fabsy1y2 > fabsy2y3) ?
    m1 * (xc - mx1) + my1 :
    m2 * (xc - mx2) + my2;
  }

  dx = x2 - xc;
  dy = y2 - yc;
  return {i: i, j: j, k: k, x: xc, y: yc, r: dx * dx + dy * dy};
}

function dedup(edges) {
  var i, j, a, b, m, n;

  for (j = edges.length; j;) {
    b = edges[--j];
    a = edges[--j];

    for (i = j; i;) {
      n = edges[--i];
      m = edges[--i];

      if ((a === m && b === n) || (a === n && b === m)) {
        edges.splice(j, 2);
        edges.splice(i, 2);
        break;
      }
    }
  }
}

function triangulate(vertices, key) {
  var n = vertices.length,
      i, j, indices, st, open, closed, edges, dx, dy, a, b, c;

  /* Bail if there aren't enough vertices to form any triangles. */
  if (n < 3)
    return [];

  /* Slice out the actual vertices from the passed objects. (Duplicate the
   * array even if we don't, though, since we need to make a supertriangle
   * later on!) */
  vertices = vertices.slice(0);
  // vertices2 = vertices.slice(0);

  if (key)
    for (i = n; i--;)
      vertices[i] = vertices[i][key];

  /* Make an array of indices into the vertex array, sorted by the
   * vertices' x-position. */
  indices = new Array(n);

  for (i = n; i--;)
    indices[i] = i;

  indices.sort(function (i, j) {
    return vertices[j][0] - vertices[i][0];
  });

  /* Next, find the vertices of the supertriangle (which contains all other
   * triangles), and append them onto the end of a (copy of) the vertex
   * array. */
  st = supertriangle(vertices);
  vertices.push(st[0], st[1], st[2]);

  /* Initialize the open list (containing the supertriangle and nothing
   * else) and the closed list (which is empty since we havn't processed
   * any triangles yet). */
  open = [circumcircle(vertices, n + 0, n + 1, n + 2)];
  closed = [];
  edges = [];

  /* Incrementally add each vertex to the mesh. */
  for (i = indices.length; i--; edges.length = 0) {
    c = indices[i];

    /* For each open triangle, check to see if the current point is
     * inside it's circumcircle. If it is, remove the triangle and add
     * it's edges to an edge list. */
    for (j = open.length; j--;) {
      /* If this point is to the right of this triangle's circumcircle,
       * then this triangle should never get checked again. Remove it
       * from the open list, add it to the closed list, and skip. */
      dx = vertices[c][0] - open[j].x;
      if (dx > 0.0 && dx * dx > open[j].r) {
        closed.push(open[j]);
        open.splice(j, 1);
        continue;
      }

      /* If we're outside the circumcircle, skip this triangle. */
      dy = vertices[c][1] - open[j].y;
      if (dx * dx + dy * dy - open[j].r > epsilon)
        continue;

      /* Remove the triangle and add it's edges to the edge list. */
      edges.push(
          open[j].i, open[j].j,
          open[j].j, open[j].k,
          open[j].k, open[j].i
      );
      open.splice(j, 1);
    }

    /* Remove any doubled edges. */
    dedup(edges);

    /* Add a new triangle for each edge. */
    for (j = edges.length; j;) {
      b = edges[--j];
      a = edges[--j];
      open.push(circumcircle(vertices, a, b, c));
    }
  }

  /* Copy any remaining open triangles to the closed list, and then
   * remove any triangles that share a vertex with the supertriangle,
   * building a list of triplets that represent triangles. */
  for (i = open.length; i--;)
    closed.push(open[i]);
  open.length = 0;

  for (i = closed.length; i--;)
    if (closed[i].i < n && closed[i].j < n && closed[i].k < n)
      open.push(closed[i].i, closed[i].j, closed[i].k);

  /* Yay, we're done! */
  var triangles = [];
  for (var ind = 0; ind < open.length; ind += 3) {
    triangles.push(
        Array.from(
            open.slice(ind, ind + 3).map(function (open_ele) {
              return vertices[open_ele];
            })
        )
    );
  }
  return triangles;
}

Diagram.prototype = {
  constructor: Diagram,

  polygons: function() {
    var edges = this.edges;

    return this.cells.map(function (cell) {
      var polygon = cell.halfedges.map(function (i) {
        return cellHalfedgeStart(cell, edges[i]);
      });
      polygon.data = cell.site.data;
      return polygon;
    });
  },

  triangles: function() {
    if (this.delaunayTriangles == null) {
      var sites = [];
      this.cells.forEach(function (cell, i) {
        sites.push(cell.site.data);
      });
      this.delaunayTriangles = triangulate(sites);
    }
    return this.delaunayTriangles;
  },

  links: function() {
    return this.edges.filter(function(edge) {
      return edge.right;
    }).map(function(edge) {
      return {
        source: edge.left.data,
        target: edge.right.data
      };
    });
  }
}
