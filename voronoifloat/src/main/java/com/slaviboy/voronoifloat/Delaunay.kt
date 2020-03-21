package com.slaviboy.voronoifloat

import android.graphics.Path
import android.graphics.PointF
import kotlin.collections.ArrayList

/**
 * Class for computing the Voronoi diagram of a set of two-dimensional points. It is based on Delaunator,
 * a fast library for computing the Delaunay triangulation using sweep algorithms. The Voronoi diagram is
 * constructed by connecting the circumcenters of adjacent triangles in the Delaunay triangulation.
 */
class Delaunay(var pointsList: ArrayList<PointF>) {

    private var _delaunator = Delaunator(pointsList)
    private var _hullIndex = IntArray(pointsList.size)
    var inedges = IntArray(pointsList.size)
    var points = _delaunator.pointsList

    lateinit var collinear: ArrayList<Int>
    lateinit var halfEdges: IntArray
    lateinit var hull: IntArray
    lateinit var triangles: IntArray

    companion object {

        fun from(coords: ArrayList<Float>): Delaunay {
            val points = ArrayList<PointF>()
            for (i in 0 until coords.size / 2) {
                points.add(PointF(coords[2 * i], coords[2 * i + 1]))
            }
            return Delaunay(points)
        }
    }

    init {
        _init()
    }

    fun _init() {
        val d = _delaunator

        // check for collinear
        if (d.hull.size > 2 && collinear(d)) {
            collinear = ArrayList<Int>()
            for (i in points.indices) {
                collinear.add(i)
            }
            getSorterIndices(points, collinear)

            val e = collinear[0]
            val f = collinear[collinear.size - 1]
            val bounds = arrayListOf<Float>(
                points[e].x,
                points[e].y,
                points[f].x,
                points[f].y
            )
            val r =
                (1e-8 * Math.sqrt(
                    0.0 +
                            (bounds[3] - bounds[1]) * (bounds[3] - bounds[1]) +
                            (bounds[2] - bounds[0]) * (bounds[2] - bounds[0])
                )).toFloat()
            val n = points.size / 2
            for (i in 0 until n) {
                points[i] = jitter(points[i], r)
            }
            this._delaunator = Delaunator(points)

        } else {
            collinear = ArrayList()
        }

        halfEdges = _delaunator.halfEdges
        hull = _delaunator.hull
        triangles = _delaunator.triangles
        inedges.fill(-1)
        _hullIndex.fill(-1)
        val inedges = inedges
        val hullIndex = _hullIndex

        // Compute an index from each point to an (arbitrary) incoming halfedge
        // Used to give the first neighbor of each point; for this reason,
        // on the hull we give priority to exterior halfEdges
        var n = halfEdges.size
        for (e in 0 until n) {
            val index = if (e % 3 == 2) {
                e - 2
            } else {
                e + 1
            }
            val p = triangles[index]
            if (halfEdges[e] == -1 || inedges[p] == -1) {
                inedges[p] = e
            }
        }

        n = hull.size
        for (i in 0 until n) {
            hullIndex[hull[i]] = i
        }

        // degenerate case: 1 or 2 (distinct) points
        if (hull.size <= 2 && hull.size > 0) {
            this.triangles = intArrayOf(-1, -1, -1)
            this.halfEdges = intArrayOf(-1, -1, -1)
            this.triangles[0] = hull[0]
            this.triangles[1] = hull[1]
            this.triangles[2] = hull[1]
            inedges[hull[0]] = 1
            if (hull.size == 2) {
                inedges[hull[1]] = 0
            }
        }
    }


    fun neighbors(i: Int) = sequence {

        // degenerate case with several collinear points
        if (collinear.size > 0) {
            val l = collinear.indexOf(i)
            if (l > 0) yield(collinear[l - 1])
            if (l < collinear.size - 1) yield(collinear[l + 1])

        } else {

            val e0 = inedges[i]
            // coincident point
            if (e0 != -1) {
                var e = e0
                var p0 = -1
                do {
                    p0 = triangles[e]
                    yield(p0)
                    e = if (e % 3 == 2) {
                        e - 2
                    } else {
                        e + 1
                    }
                    // bad triangulation
                    if (triangles[e] != i) {
                        break
                    }
                    e = halfEdges[e]
                    if (e == -1) {
                        val p = hull[(_hullIndex[i] + 1) % hull.size]
                        if (p != p0) {
                            yield(p)
                        }
                        break
                    }
                } while (e != e0)
            }
        }
    }

    fun find(p: PointF, i: Int = 0): Int {
        if (p.x.isNaN() || p.y.isNaN()) {
            return -1
        }
        var _i = i
        val i0 = _i
        var c = _step(p, _i)
        while (c >= 0 && c != _i && c != i0) {
            c = _step(p, _i)
            _i = c
        }
        return c
    }

    fun _step(p: PointF, i: Int): Int {
        if (inedges[i] == -1 || points.size == 0) {
            return (i + 1) % (points.size shr 1)
        }
        var c = i;
        var dc =
            (p.x - points[i].x) * (p.x - points[i].x) +
                    (p.y - points[i].y) * (p.y - points[i].y)
        val e0 = inedges[i]
        var e = e0
        do {
            val t = triangles[e]
            val dt = (p.x - points[t].x) * (p.x - points[t].x) +
                    (p.y - points[t].y) * (p.y - points[t].y)
            if (dt < dc) {
                dc = dt
                c = t
            }
            e = if (e % 3 == 2) {
                e - 2
            } else {
                e + 1
            }
            // bad triangulation
            if (triangles[e] != i) {
                break
            }
            e = halfEdges[e]
            if (e == -1) {
                e = hull[(_hullIndex[i] + 1) % hull.size]
                if (e != t) {
                    if ((p.x - points[e].x) * (p.x - points[e].x) +
                        (p.y - points[e].y) * (p.y - points[e].y) < dc
                    ) {
                        return e
                    }
                }
                break
            }
        } while (e != e0)
        return c
    }

    fun render(path: Path): Path {
        val n = halfEdges.size
        for (i in 0 until n) {
            val j = halfEdges[i]
            if (j < i) {
                continue
            }
            val ti = triangles[i]
            val tj = triangles[j]
            path.moveTo(points[ti].x.toFloat(), points[ti].y.toFloat())
            path.lineTo(points[tj].x.toFloat(), points[tj].y.toFloat())

        }
        renderHull(path)
        return path
    }

    fun renderHull(path: Path): Path {
        val h = hull[0]
        val n = hull.size
        path.moveTo(points[h].x.toFloat(), points[h].y.toFloat())
        for (i in 0 until n) {
            val h = hull[i]
            path.lineTo(points[h].x.toFloat(), points[h].y.toFloat())
        }
        path.close()
        return path
    }

    fun renderPoints(path: Path, r: Int = 3, range: IntRange = (0 until points.size)): Path {

        for (i in range) {
            val x = points[i].x
            val y = points[i].y
            path.moveTo((x + r).toFloat(), y.toFloat())
            path.addCircle(
                x.toFloat(), y.toFloat(), r.toFloat(), Path.Direction.CW
            )
        }
        return path
    }

    fun renderTriangle(i: Int, path: Path): Path {
        val _i = i * 3
        val t0 = triangles[_i]
        val t1 = triangles[_i + 1]
        val t2 = triangles[_i + 2]
        path.moveTo(points[t0].x.toFloat(), points[t0].y.toFloat())
        path.lineTo(points[t1].x.toFloat(), points[t1].y.toFloat())
        path.lineTo(points[t2].x.toFloat(), points[t2].y.toFloat())
        path.close()
        return path
    }

    fun triangleCenter(i: Int): PointF {
        val _i = i * 3
        val t0 = triangles[_i]
        val t1 = triangles[_i + 1]
        val t2 = triangles[_i + 2]

        return PointF(
            (points[t0].x + points[t1].x + points[t2].x) / 3,
            (points[t0].y + points[t1].y + points[t2].y) / 3
        )
    }

    fun trianglePolygons() = sequence {
        for (i in 0 until triangles.size / 3) {
            yield(trianglePolygon(i))
        }
    }

    fun trianglePolygon(i: Int): Path {
        val path = Path()
        this.renderTriangle(i, path)
        return path
    }

    fun update(): Delaunay {
        _delaunator.update()
        _init()
        return this
    }

    fun getSorterIndices(a: ArrayList<PointF>, b: ArrayList<Int>) {
        val temp = a.copy()
        for (x in a.indices) {
            for (y in x + 1 until a.size - 1) {
                if (temp[y].x < temp[x].x || (temp[y].x == temp[x].x && temp[y].y < temp[x].y)) {
                    temp.swap(y, x)
                    b.swap(x, y)
                }
            }
        }
    }

    fun <T> ArrayList<T>.copy(): ArrayList<T> {
        val temp = ArrayList<T>()
        for (i in 0 until this.size) {
            temp.add(this[i])
        }
        return temp
    }

    fun <T> ArrayList<T>.swap(index1: Int, index2: Int) {
        val tmp = this[index1] // 'this' corresponds to the list
        this[index1] = this[index2]
        this[index2] = tmp
    }

    // A triangulation is collinear if all its triangles have a non-null area
    fun collinear(d: Delaunator): Boolean {

        val triangles = d.triangles
        val coords = d.pointsList
        for (i in triangles.indices step 3) {
            val a = triangles[i]
            val b = triangles[i + 1]
            val c = triangles[i + 2]
            val cross = (coords[c].x - coords[a].x) * (coords[b].y - coords[a].y)
            -(coords[b].x - coords[a].x) * (coords[c].y - coords[a].y)
            if (cross > 1e-10) {
                return false
            }
        }
        return true
    }

    fun jitter(p: PointF, r: Float): PointF {
        val x = (p.x + Math.sin(0.0 + p.x + p.y) * r).toFloat()
        val y = (p.y + Math.cos(0.0 + p.x - p.y) * r).toFloat()
        return PointF(x, y)
    }
}