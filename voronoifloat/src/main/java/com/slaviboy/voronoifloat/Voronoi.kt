package com.slaviboy.voronoifloat

import android.graphics.Path
import android.graphics.PointF
import android.graphics.RectF
import android.util.Log

class Voronoi(var delaunay: Delaunay, var bound: RectF) {

    private lateinit var vectors: Array<PointF>
      lateinit var circumcenters: Array<PointF>

    init {
        if (bound.right.isNaN() || bound.bottom.isNaN() || bound.right < bound.left || bound.bottom < bound.top) {
            throw IllegalArgumentException("Invalid bounds!")
        }
        _init()
    }

    private fun _init() {

        val points = delaunay.points
        val hull = delaunay.hull
        val triangles = delaunay.triangles

        // compute circumcenters
        circumcenters = Array<PointF>(triangles.size / 3) { PointF() }
        for (i in triangles.indices step 3) {

            val t1 = triangles[i]
            val p1 = points[t1]

            val t2 = triangles[i + 1]
            val p2 = points[t2]

            val t3 = triangles[i + 2]
            val p3 = points[t3]

            val dx = p2.x - p1.x
            val dy = p2.y - p1.y
            val ex = p3.x - p1.x
            val ey = p3.y - p1.y
            val bl = dx * dx + dy * dy
            val cl = ex * ex + ey * ey
            val ab = (dx * ey - dy * ex) * 2


            circumcenters[i / 3] = if (ab.isNaN() || ab == 0.0f) {
                // degenerate case (collinear diagram)
                val x = ((p1.x + p3.x) / 2 - 1e8 * ey).toFloat()
                val y = ((p1.y + p3.y) / 2 + 1e8 * ex).toFloat()
                PointF(x, y)
            } else if (Math.abs(ab) < 1e-8) {
                // almost equal points (degenerate triangle)
                PointF(
                    (p1.x + p3.x) / 2f,
                    (p1.y + p3.y) / 2f
                )
            } else {
                val d = 1 / ab
                PointF(
                    p1.x + (ey * bl - dy * cl) * d,
                    p1.y + (dx * cl - ex * bl) * d
                )
            }
        }


        vectors = Array<PointF>(delaunay.points.size * 2) { PointF() }

        // Compute exterior cell rays.
        var h = hull[hull.size - 1]
        var i0: Int
        var i1 = h * 2
        var p0: PointF
        var p1 = points[h]
        for (i in hull.indices) {
            h = hull[i]
            i0 = i1
            i1 = h * 2
            p0 = PointF(p1.x, p1.y)
            p1 = points[h]

            vectors[i0 + 1] = PointF(p0.y - p1.y, p1.x - p0.x)
            vectors[i1] = PointF(p0.y - p1.y, p1.x - p0.x)
        }
    }

    fun update(): Voronoi {
        delaunay.update()
        _init()
        return this
    }

    fun render(path: Path): Path {
        val halfedges = delaunay.halfEdges
        val inedges = delaunay.inedges
        val hull = delaunay.hull

        if (hull.size <= 1) {
            return path
        }

        for (i in halfedges.indices) {
            val j = halfedges[i]
            if (j < i) {
                continue
            }
            renderSegment(circumcenters[i / 3], circumcenters[j / 3], path)
        }

        var h0: Int
        var h1 = hull[hull.size - 1]
        for (i in hull.indices) {
            h0 = h1
            h1 = hull[i]

            val t = inedges[h1] / 3
            val v = h0 * 2
            val p = project(circumcenters[t], vectors[v + 1])
            if (p != null) {
                renderSegment(circumcenters[t], p, path)
            }
        }
        return path
    }

    fun renderBounds(path: Path): Path {
        path.addRect(
            bound.left.toFloat(),
            bound.top.toFloat(),
            bound.right.toFloat(),
            bound.bottom.toFloat(),
            Path.Direction.CW
        )
        return path
    }

    fun renderCell(i: Int, path: Path): Path {
        val points = clip(i)

        if (points != null) {
            path.moveTo(points[0].x.toFloat(), points[0].y.toFloat())
            var n = points.size
            while (points[0].x == points[n - 1].x && points[1].y == points[n - 1].y && n > 1) {
                n -= 1
            }

            var i = 1
            while (i < n) {
                if (points[i].x != points[i - 1].x || points[i].y != points[i - 1].y) {
                    path.lineTo(points[i].x.toFloat(), points[i].y.toFloat())
                }
                i++
            }
            path.close()
        }
        return path
    }

    fun getCell(i: Int): ArrayList<PointF> {
        return clip(i) ?: ArrayList()
    }

    /**
     * Lloyd relaxation for all cells in the current voronoi
     * @param iterations how many time to do the relaxation algorithm
     */
    fun relax(iterations: Int) {
        for (t in 0 until iterations) {
            for (i in 0 until delaunay.points.size) {

                val pointsArray = getCell(i)
                val centroid = Polygon.centroid(pointsArray)
                if (centroid.x > 0 && centroid.x < bound.width() &&
                    centroid.y > 0 && centroid.y < bound.height()
                ) {
                    delaunay.points[i] = centroid
                }
            }
            update()
        }
    }

    private fun clip(i: Int): ArrayList<PointF>? {
        // degenerate case (1 valid point: return the box)
        if (i == 0 && delaunay.hull.size == 1) {
            return arrayListOf(
                PointF(bound.right, bound.top), PointF(bound.right, bound.bottom),
                PointF(bound.left, bound.bottom), PointF(bound.left, bound.top)
            )
        }

        val points = cell(i) ?: return null
        val v = i * 2
        return if (vectors[v].x != 0.0f || vectors[v].y != 0.0f) {
            clipInfinite(i, points, vectors[v], vectors[v + 1])
        } else {
            clipFinite(i, points)
        }
    }

    private fun cell(i: Int): ArrayList<PointF> {
        val e0 = delaunay.inedges[i]
        if (e0 == -1) return ArrayList()  // coincident point

        val points = ArrayList<PointF>()
        var e = e0
        do {
            val t = e / 3
            points.add(circumcenters[t])
            e = if (e % 3 == 2) {
                e - 2
            } else {
                e + 1
            }
            if (delaunay.triangles[e] != i) break; // bad triangulation
            e = delaunay.halfEdges[e]
        } while (e != e0 && e != -1)
        return points
    }

    private fun clipInfinite(
        i: Int,
        points: ArrayList<PointF>,
        vp0: PointF,
        vpn: PointF
    ): ArrayList<PointF> {
        var pts = ArrayList<PointF>(points)
        var p: PointF? = null

        p = project(pts[0], vp0)
        if (p != null) {
            pts.add(0, p)
        }

        p = project(pts[pts.size - 1], vpn)
        if (p != null) {
            pts.add(p)
        }

        val _pts = clipFinite(i, pts)
        if (_pts != null) {
            pts = _pts

            var n = pts.size
            var c0: Int
            var c1 = edgecode(pts[n - 1])

            var j = 0
            while (j < n) {
                c0 = c1
                c1 = edgecode(pts[j])
                if (c0 != 0 && c1 != 0) {
                    j = edge(i, c0, c1, pts, j)
                    n = pts.size
                }
                j++
            }
        } else {
            if (contains(
                    i,
                    PointF((bound.left + bound.right) / 2, (bound.top + bound.bottom) / 2)
                )
            ) {
                pts = arrayListOf(
                    PointF(bound.left, bound.top),
                    PointF(bound.right, bound.top),
                    PointF(bound.right, bound.bottom),
                    PointF(bound.left, bound.bottom)
                )
            }
        }
        return pts
    }

    private fun clipFinite(i: Int, points: ArrayList<PointF>): ArrayList<PointF>? {

        var pts: ArrayList<PointF>? = null
        var p0: PointF
        var p1 = points[points.size - 1]
        var c0: Int
        var c1 = regioncode(p1)
        var e0: Int
        var e1 = 0

        for (j in 0 until points.size) {
            p0 = PointF(p1.x, p1.y)
            p1 = points[j]
            c0 = c1
            c1 = regioncode(p1)
            if (c0 == 0 && c1 == 0) {
                e0 = e1
                e1 = 0
                if (pts != null) {
                    pts.add(p1)
                } else {
                    pts = arrayListOf(p1)
                }
            } else {
                var s: Array<PointF>?
                var sp0: PointF
                var sp1: PointF
                if (c0 == 0) {
                    s = clipSegment(p0, p1, c0, c1)
                    if (s == null) continue
                    sp0 = s[0]
                    sp1 = s[1]
                } else {
                    s = clipSegment(p1, p0, c1, c0)
                    if (s == null) continue
                    sp1 = s[0]
                    sp0 = s[1]
                    e0 = e1
                    e1 = edgecode(sp0)
                    if (e0 != 0 && e1 != 0) {
                        if (pts != null) {
                            edge(i, e0, e1, pts, pts.size)
                        }
                    }
                    if (pts != null) {
                        pts.add(sp0)
                    } else {
                        pts = arrayListOf(sp0)
                    }
                }
                e0 = e1
                e1 = edgecode(sp1)
                if (e0 != 0 && e1 != 0) {
                    if (pts != null) {
                        edge(i, e0, e1, pts, pts.size)
                    }
                }
                if (pts != null) {
                    pts.add(sp1)
                } else {
                    pts = arrayListOf(sp1)
                }
            }
        }

        if (pts != null) {
            e0 = e1
            e1 = edgecode(pts[0])
            if (e0 != 0 && e1 != 0) edge(i, e0, e1, pts, pts.size)
        } else if (contains(
                i,
                PointF((bound.left + bound.right) / 2, (bound.top + bound.bottom) / 2)
            )
        ) {
            return arrayListOf(
                PointF(bound.right, bound.top),
                PointF(bound.right, bound.bottom),
                PointF(bound.left, bound.bottom),
                PointF(bound.left, bound.top)
            )
        }
        return pts
    }

    private fun edge(i: Int, e0: Int, e1: Int, pts: ArrayList<PointF>, j: Int): Int {
        var _j = j
        var _e0 = e0
        loop@ while (_e0 != e1) {
            lateinit var p: PointF
            when (_e0) {
                0b0101 -> {
                    // top-left
                    _e0 = 0b0100
                    continue@loop
                }
                0b0100 -> {
                    // top
                    _e0 = 0b0110
                    p = PointF(bound.right, bound.top)
                }
                0b0110 -> {
                    // top-right
                    _e0 = 0b0010
                    continue@loop
                }
                0b0010 -> {
                    // right
                    _e0 = 0b1010
                    p = PointF(bound.right, bound.bottom)
                }
                0b1010 -> {
                    // bottom-right
                    _e0 = 0b1000
                    continue@loop
                }
                0b1000 -> {
                    // bottom
                    _e0 = 0b1001
                    p = PointF(bound.left, bound.bottom)
                }
                0b1001 -> {
                    // bottom-left
                    _e0 = 0b0001
                    continue@loop
                }
                0b0001 -> {
                    // left
                    _e0 = 0b0101
                    p = PointF(bound.left, bound.top)
                }
            }
            if ((pts.size <= _j || pts[_j].x != p.x || pts[_j].y != p.y) && contains(i, p)) {
                pts.add(_j, p)
                _j++
            }
        }

        if (pts.size > 2) {
            var i = 0
            while (i < pts.size) {
                val j = (i + 1) % pts.size
                val k = (i + 2) % pts.size

                if ((pts[i].x == pts[j].x && pts[j].x == pts[k].x) || (pts[i].y == pts[j].y && pts[j].y == pts[k].y)) {
                    for (t in j until (j + 1)) {
                        pts.removeAt(j)
                    }
                    i -= 1
                }
                i += 1
            }
        }
        return _j
    }

    private fun renderSegment(p0: PointF, p1: PointF, path: Path) {
        val c0 = regioncode(p0)
        val c1 = regioncode(p1)
        if (c0 == 0 && c1 == 0) {
            path.moveTo(p0.x.toFloat(), p0.y.toFloat())
            path.lineTo(p1.x.toFloat(), p1.y.toFloat())
        } else {
            val s = clipSegment(p0, p1, c0, c1)
            if (s != null) {
                path.moveTo(s[0].x.toFloat(), s[0].y.toFloat())
                path.lineTo(s[1].x.toFloat(), s[1].y.toFloat())
            }
        }
    }

    private fun project(p0: PointF, vp: PointF): PointF? {
        var t = Float.POSITIVE_INFINITY
        var p = PointF()
        var c: Float

        if (vp.y < 0) {
            // top
            if (p0.y <= bound.top) return null
            c = (bound.top - p0.y) / vp.y
            if (c < t) {
                t = c
                p = PointF(p0.x + t * vp.x, bound.top)
            }
        } else if (vp.y > 0) {
            // bottom
            if (p0.y >= bound.bottom) return null
            c = (bound.bottom - p0.y) / vp.y
            if (c < t) {
                t = c
                p = PointF(p0.x + t * vp.x, bound.bottom)
            }
        }

        if (vp.x > 0) {
            // right
            if (p0.x >= bound.right) return null
            c = (bound.right - p0.x) / vp.x
            if (c < t) {
                t = c
                p = PointF(bound.right, p0.y + t * vp.y)
            }
        } else if (vp.x < 0) {
            // left
            if (p0.x <= bound.left) return null
            c = (bound.left - p0.x) / vp.x
            if (c < t) {
                t = c
                p = PointF(bound.left, p0.y + t * vp.y)
            }
        }

        return p
    }

    private fun contains(i: Int, p: PointF): Boolean {
        if (p.x.isNaN() || p.y.isNaN()) return false
        return delaunay._step(p, i) == i
    }

    private fun edgecode(p: PointF): Int {
        val b1 = if (p.x == bound.left) {
            0b0001
        } else {
            if (p.x == bound.right) {
                0b0010
            } else {
                0b0000
            }
        }

        val b2 = if (p.y == bound.top) {
            0b0100
        } else {
            if (p.y == bound.bottom) {
                0b1000
            } else {
                0b0000
            }
        }
        return b1 or b2
    }

    private fun regioncode(p: PointF): Int {
        val b1 = if (p.x < bound.left) {
            0b0001
        } else {
            if (p.x > bound.right) 0b0010
            else 0b0000
        }

        val b2 = if (p.y < bound.top) {
            0b0100
        } else {
            if (p.y > bound.bottom) 0b1000
            else 0b0000
        }
        return b1 or b2
    }

    private fun clipSegment(p0: PointF, p1: PointF, c0: Int, c1: Int): Array<PointF>? {

        var _c0 = c0
        var _c1 = c1
        var _p0 = p0
        var _p1 = p1
        while (true) {
            if (_c0 == 0 && _c1 == 0) {
                return arrayOf(_p0, _p1)
            }

            if ((_c0 and _c1) != 0) return null
            val c = if (_c0 != 0) {
                _c0
            } else {
                _c1
            }

            val p = when {
                (c and 0b1000) != 0 -> {
                    PointF(
                        _p0.x + (_p1.x - _p0.x) * (bound.bottom - _p0.y) / (_p1.y - _p0.y),
                        bound.bottom
                    )
                }
                (c and 0b0100) != 0 -> {
                    PointF(
                        _p0.x + (_p1.x - _p0.x) * (bound.top - _p0.y) / (_p1.y - _p0.y),
                        bound.top
                    )
                }
                (c and 0b0010) != 0 -> {
                    PointF(
                        bound.right,
                        _p0.y + (_p1.y - _p0.y) * (bound.right - _p0.x) / (_p1.x - _p0.x)
                    )
                }
                else -> {
                    PointF(
                        bound.left,
                        _p0.y + (_p1.y - _p0.y) * (bound.left - _p0.x) / (_p1.x - _p0.x)
                    )
                }
            }

            if (_c0 != 0) {
                _p0 = p
                _c0 = regioncode(p)
            } else {
                _p1 = p
                _c1 = regioncode(p)
            }
        }
    }
}