package com.slaviboy.voronoifloat

import android.graphics.PointF
import android.util.Log
import kotlin.math.*

/**
 *  Class that construct the Voronoi diagram by getting the Delaunay triangulation
 *  with all the circumcircles and their centers of the same set of points. Then connecting
 *  the centers of the circumcircles to produce the Voronoi diagram. Learn more
 *  https://en.wikipedia.org/wiki/Delaunay_triangulation
 *
 *  Inspired by the JS implementation
 *  https://github.com/mapbox/delaunator
 *
 *  @param pointsList array with points
 */
class Delaunator(var pointsList: ArrayList<PointF>) {

    lateinit var triangles: IntArray
    lateinit var halfEdges: IntArray
    lateinit var hull: IntArray

    private var _triangles: IntArray
    private var _halfEdges: IntArray
    private var _hashSize: Int
    private var _hullPrev: IntArray
    private var _hullNext: IntArray
    private var _hullTri: IntArray
    private var _hullHash: IntArray
    private var _ids: IntArray
    private var _dists: FloatArray
    private var _cx: Float = 0.0f
    private var _cy: Float = 0.0f
    private var _hullStart: Int = 0
    private var trianglesLen: Int = 0

    // static methods and variable
    companion object {
        val EPSILON = Math.pow(2.0, -52.0).toFloat()
        val EDGE_STACK = IntArray(512)

        /**
         * Init delaunator object from array with coordinates
         * arrayListOf(x1,y1, x2,y2, x3,y3, ....)
         *
         * @param coordinates array with coordinates
         */
        fun from(coordinates: ArrayList<Float>): Delaunator {

            // convert to array with points
            val points = ArrayList<PointF>()
            for (i in 0 until coordinates.size / 2) {
                points.add(PointF(coordinates[2 * i], coordinates[2 * i + 1]))
            }
            return Delaunator(points)
        }

    }

    init {
        val n = pointsList.size

        // arrays that will store the triangulation graph
        val maxTriangles = max(2 * n - 5, 0)
        _triangles = IntArray(maxTriangles * 3)
        _halfEdges = IntArray(maxTriangles * 3)

        // temporary arrays for tracking the edges of the advancing convex hull
        _hashSize = ceil(sqrt(n.toFloat())).toInt()
        _hullPrev = IntArray(n)                     // edge to prev edge
        _hullNext = IntArray(n)                     // edge to next edge
        _hullTri = IntArray(n)                      // edge to adjacent triangle
        _hullHash = IntArray(_hashSize)             // angular edge hash
        _hullHash.fill(-1, 0, _hashSize)

        // temporary arrays for sorting points
        _ids = IntArray(n)
        _dists = FloatArray(n)

        update()
    }

    fun update() {

        // populate an array of point indices  calculate input data bbox
        var minX = Float.POSITIVE_INFINITY
        var minY = Float.POSITIVE_INFINITY
        var maxX = Float.NEGATIVE_INFINITY
        var maxY = Float.NEGATIVE_INFINITY

        for (i: Int in pointsList.indices) {
            val x = pointsList[i].x
            val y = pointsList[i].y
            if (x < minX) minX = x
            if (y < minY) minY = y
            if (x > maxX) maxX = x
            if (y > maxY) maxY = y
            this._ids[i] = i
        }
        val c = PointF((minX + maxX) / 2.0f, (minY + maxY) / 2.0f)

        var minDist = Float.POSITIVE_INFINITY
        var i0 = 0
        var i1 = 0
        var i2 = 0

        // pick a seed point close to the center
        for (i: Int in pointsList.indices) {
            val d: Float = dist(c, pointsList[i])
            if (d < minDist) {
                i0 = i
                minDist = d
            }
        }
        val p0 = PointF(pointsList[i0].x, pointsList[i0].y)
        minDist = Float.POSITIVE_INFINITY

        // find the point closest to the seed
        for (i: Int in pointsList.indices) {
            if (i == i0) {
                continue
            }
            val d = dist(p0, pointsList[i])
            if (d < minDist && d > 0) {
                i1 = i
                minDist = d
            }
        }

        val p1 = PointF(pointsList[i1].x, pointsList[i1].y)
        var minRadius = Float.POSITIVE_INFINITY

        // find the third point which forms the smallest circumcircle with the first two
        for (i: Int in pointsList.indices) {
            if (i == i0 || i == i1) {
                continue
            }
            val r: Float = circumradius(p0, p1, pointsList[i])
            if (r < minRadius) {
                i2 = i
                minRadius = r
            }
        }

        val p2 = PointF(pointsList[i2].x, pointsList[i2].y)

        if (minRadius == Float.POSITIVE_INFINITY) {
            // order collinear points by dx (or dy if all x are identical)
            // and return the list as a hull
            for (i: Int in pointsList.indices) {
                _dists[i] =
                    if (pointsList[i].x - pointsList[0].x > 0) {
                        pointsList[i].x - pointsList[0].x
                    } else {
                        pointsList[i].y - pointsList[1].y
                    }
            }

            quicksort(_ids, _dists, 0, pointsList.size - 1)

            val hullTemp = IntArray(pointsList.size)
            var j = 0
            var d0 = Float.NEGATIVE_INFINITY
            for (i: Int in pointsList.indices) {
                val id = _ids[i]
                if (_dists[id] > d0) {
                    hullTemp[j++] = _ids[i]
                    d0 = _dists[id]
                }
            }

            hull = hullTemp.copyOfRange(0, j)
            triangles = IntArray(0)
            halfEdges = IntArray(0)
            return
        }

        // swap the order of the seed points for counter-clockwise orientation
        if (orient(p0, p1, p2)) {
            val i = i1
            val x = p1.x
            val y = p1.y

            i1 = i2
            p1.x = p2.x
            p1.y = p2.y

            i2 = i
            p2.x = x
            p2.y = y
        }

        val center = circumcenter(p0, p1, p2)
        _cx = center.x
        _cy = center.y

        for (i: Int in pointsList.indices) {
            _dists[i] = dist(pointsList[i], center)
        }

        // sort the points by distance from the seed triangle circumcenter
        quicksort(_ids, _dists, 0, pointsList.size - 1)

        // set up the seed triangle as the starting hull
        _hullStart = i0
        var hullSize = 3

        _hullNext[i0] = i1
        _hullPrev[i2] = i1
        _hullNext[i1] = i2
        _hullPrev[i0] = i2
        _hullNext[i2] = i0
        _hullPrev[i1] = i0

        _hullTri[i0] = 0
        _hullTri[i1] = 1
        _hullTri[i2] = 2

        _hullHash.fill(-1)
        _hullHash[hashKey(p0).toInt()] = i0
        _hullHash[hashKey(p1).toInt()] = i1
        _hullHash[hashKey(p2).toInt()] = i2

        trianglesLen = 0
        addTriangle(i0, i1, i2, -1, -1, -1)

        var xp = 0.0f
        var yp = 0.0f
        for (k in _ids.indices) {
            val i = _ids[k]
            val x = pointsList[i].x
            val y = pointsList[i].y

            // skip near-duplicate points
            if (k > 0 && abs(x - xp) <= EPSILON && abs(y - yp) <= EPSILON) {
                continue
            }
            xp = x
            yp = y

            // skip seed triangle points
            if (i == i0 || i == i1 || i == i2) {
                continue
            }

            // find a visible edge on the convex hull using edge hash
            var start = 0
            val key = hashKey(pointsList[i]).toInt()
            for (j in 0 until _hashSize) {
                start = _hullHash[(key + j) % _hashSize]
                if (start != -1 && start != _hullNext[start]) {
                    break
                }
            }

            start = _hullPrev[start]
            var e = start
            var q = _hullNext[e]
            while (!orient(pointsList[i], pointsList[e], pointsList[q])) {
                e = q
                if (e == start) {
                    e = -1
                    break
                }
                q = _hullNext[e]
            }

            // likely a near-duplicate point- skip it
            if (e == -1) {
                continue
            }

            // add the first triangle from the point
            var t = addTriangle(e, i, _hullNext[e], -1, -1, _hullTri[e])

            // recursively flip triangles from the point until they satisfy the Delaunay condition
            _hullTri[i] = legalize(t + 2)
            _hullTri[e] = t // keep track of boundary triangles on the hull
            hullSize++

            // walk forward through the hull, adding more triangles and flipping recursively
            var n = _hullNext[e]
            q = _hullNext[n]
            while (orient(pointsList[i], pointsList[n], pointsList[q])) {
                t = addTriangle(n, i, q, _hullTri[i], -1, _hullTri[n])
                _hullTri[i] = legalize(t + 2)
                _hullNext[n] = n // mark as removed
                hullSize--
                n = q
                q = _hullNext[n]
            }

            // walk backward from the other side, adding more triangles and flipping
            if (e == start) {
                q = _hullPrev[e]
                while (orient(pointsList[i], pointsList[q], pointsList[e])) {
                    t = addTriangle(q, i, e, -1, _hullTri[e], _hullTri[q])
                    legalize(t + 2)
                    _hullTri[q] = t
                    _hullNext[e] = e // mark as removed
                    hullSize--
                    e = q
                    q = _hullPrev[e]
                }
            }

            // update the hull indices
            _hullStart = e
            _hullPrev[i] = e
            _hullNext[e] = i
            _hullPrev[n] = i
            _hullNext[i] = n

            // save the two new edges in the hash table
            _hullHash[hashKey(pointsList[i]).toInt()] = i
            _hullHash[hashKey(pointsList[e]).toInt()] = e
        }

        hull = IntArray(hullSize)
        var e = _hullStart
        for (i in 0 until hullSize) {
            hull[i] = e
            e = _hullNext[e]
        }

        // trim typed triangle mesh arrays
        triangles = _triangles.copyOfRange(0, trianglesLen)
        halfEdges = _halfEdges.copyOfRange(0, trianglesLen)
    }

    /**
     * Monotonically increases with real angle, but doesn't need expensive trigonometry
     */
    private fun pseudoAngle(dx: Float, dy: Float): Float {
        val v = dx / (abs(dx) + abs(dy))

        // [0..1]
        return if (dy > 0f) {
            (3 - v) / 4.0f
        } else {
            (1 + v) / 4.0f
        }
    }

    /**
     * Represent the distance between two points, without the square root of both sides
     *
     * @param a first point
     * @param b second point
     * @return distance between the two points
     */
    private fun dist(a: PointF, b: PointF): Float {
        val dx = a.x - b.x
        val dy = a.y - b.y
        return dx * dx + dy * dy
    }

    /**
     * Find the circumradius, from the three points of a triangle
     */
    private fun circumradius(a: PointF, b: PointF, c: PointF): Float {
        val dx = b.x - a.x
        val dy = b.y - a.y
        val ex = c.x - a.x
        val ey = c.y - a.y

        val bl = dx * dx + dy * dy
        val cl = ex * ex + ey * ey
        val d = 0.5f / (dx * ey - dy * ex)

        val x = (ey * bl - dy * cl) * d
        val y = (dx * cl - ex * bl) * d

        return x * x + y * y
    }

    /**
     * Quicksort for rearranging the ids, by there corresponding distances
     * from smaller to bigger
     *
     * @param ids list with ids corresponding to each distance
     * @param dists list with distances corresponding to each id
     * @param left start position
     * @param right end position
     *
     *
     * @example
     *  --------------------
     *  input:
     *  --------------------
     *  ids = [0,1,2,3,4,5]
     *  dists = [24,5,23,12,77,52]
     *  left = 0
     *  right = ids.size
     *  --------------------
     *  =>  ids = [1, 3, 2, 0, 5, 4]
     */
    private fun quicksort(ids: IntArray, dists: FloatArray, left: Int, right: Int) {
        if (right - left <= 20) {
            for (i: Int in (left + 1)..right) {
                val temp = ids[i]
                val tempDist = dists[temp]
                var j = i - 1
                while (j >= left && dists[ids[j]] > tempDist)
                    ids[j + 1] = ids[j--]
                ids[j + 1] = temp
            }
        } else {
            val median = (left + right) shr 1
            var i = left + 1
            var j = right
            swap(ids, median, i)
            if (dists[ids[left]] > dists[ids[right]]) swap(ids, left, right)
            if (dists[ids[i]] > dists[ids[right]]) swap(ids, i, right)
            if (dists[ids[left]] > dists[ids[i]]) swap(ids, left, i)

            val temp = ids[i]
            val tempDist = dists[temp]
            while (true) {
                do i++ while (dists[ids[i]] < tempDist)
                do j-- while (dists[ids[j]] > tempDist)
                if (j < i) break
                swap(ids, i, j)
            }
            ids[left + 1] = ids[j]
            ids[j] = temp

            if (right - i + 1 >= j - left) {
                quicksort(ids, dists, i, right)
                quicksort(ids, dists, left, j - 1)
            } else {
                quicksort(ids, dists, left, j - 1)
                quicksort(ids, dists, i, right)
            }
        }
    }

    /**
     * Swap elements in array by given indices
     *
     * @param arr array with integers
     * @param i first index
     * @param j second index
     */
    private fun swap(arr: IntArray, i: Int, j: Int) {
        val tmp = arr[i]
        arr[i] = arr[j]
        arr[j] = tmp
    }

    // a more robust orientation test that's stable in a given triangle (to fix robustness issues)
    private fun orient(r: PointF, q: PointF, p: PointF): Boolean {
        val sign = when {
            orientIfSure(p, r, q) > 0 -> {
                orientIfSure(p, r, q)
            }
            orientIfSure(r, q, p) > 0 -> {
                orientIfSure(r, q, p)
            }
            else -> {
                orientIfSure(q, p, r)
            }
        }
        return sign < 0
    }

    /**
     * Return 2d orientation sign if we're confident in it through J. Shewchuk's error bound check
     */
    private fun orientIfSure(p: PointF, r: PointF, q: PointF): Float {
        val l = (r.y - p.y) * (q.x - p.x)
        val m = (r.x - p.x) * (q.y - p.y)
        return if (abs(l - m) >= 3.3306690738754716e-16 * abs(l + m)) {
            l - m
        } else {
            0.0f
        }
    }

    /**
     * Find the circumcircle center of a triangle with given points
     *
     * @param a first point
     * @param b second point
     * @param c third point
     */
    private fun circumcenter(a: PointF, b: PointF, c: PointF): PointF {

        val dx = b.x - a.x
        val dy = b.y - a.y
        val bl = dx * dx + dy * dy

        val ex = c.x - a.x
        val ey = c.y - a.y
        val cl = ex * ex + ey * ey

        val d = 0.5f / (dx * ey - dy * ex)
        val x = a.x + (ey * bl - dy * cl) * d
        val y = a.y + (dx * cl - ex * bl) * d

        return PointF(x, y)
    }

    private fun hashKey(p: PointF): Float {
        return floor(
            pseudoAngle(p.x - _cx, p.y - _cy) * _hashSize
        ) % _hashSize
    }

    private fun addTriangle(i0: Int, i1: Int, i2: Int, a: Int, b: Int, c: Int): Int {
        val t = trianglesLen

        _triangles[t] = i0
        _triangles[t + 1] = i1
        _triangles[t + 2] = i2

        link(t, a)
        link(t + 1, b)
        link(t + 2, c)

        trianglesLen += 3
        return t
    }

    private fun link(a: Int, b: Int) {
        _halfEdges[a] = b
        if (b != -1) _halfEdges[b] = a
    }

    private fun legalize(_a: Int): Int {

        var a = _a
        var i = 0
        var ar: Int

        // recursion eliminated with a fixed-size stack
        while (true) {
            val b = _halfEdges[a]

            /* if the pair of triangles doesn't satisfy the Delaunay condition
             * (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
             * then do the same check/flip recursively for the new pair of triangles
             *
             *           pl                    pl
             *          /||\                  /  \
             *       al/ || \bl            al/    \a
             *        /  ||  \              /      \
             *       /  a||b  \    flip    /___ar___\
             *     p0\   ||   /p1   =>   p0\---bl---/p1
             *        \  ||  /              \      /
             *       ar\ || /br             b\    /br
             *          \||/                  \  /
             *           pr                    pr
             */
            val a0 = a - a % 3
            ar = a0 + (a + 2) % 3

            // convex hull edge
            if (b == -1) {
                if (i == 0) {
                    break
                }
                a = EDGE_STACK[--i]
                continue
            }

            val b0 = b - b % 3
            val al = a0 + (a + 1) % 3
            val bl = b0 + (b + 2) % 3

            val p0 = _triangles[ar]
            val pr = _triangles[a]
            val pl = _triangles[al]
            val p1 = _triangles[bl]

            val illegal = inCircle(
                pointsList[p0], pointsList[pr],
                pointsList[pl], pointsList[p1]
            )

            if (illegal) {
                _triangles[a] = p1
                _triangles[b] = p0

                val hbl = _halfEdges[bl]

                // edge swapped on the other side of the hull (rare)  fix the halfedge reference
                if (hbl == -1) {
                    var e = _hullStart
                    do {
                        if (_hullTri[e] == bl) {
                            _hullTri[e] = a
                            break
                        }
                        e = _hullPrev[e]
                    } while (e != _hullStart)
                }
                link(a, hbl)
                link(b, _halfEdges[ar])
                link(ar, bl)

                val br = b0 + (b + 1) % 3

                // don't worry about hitting the cap: it can only happen on extremely degenerate input
                if (i < EDGE_STACK.size) {
                    EDGE_STACK[i++] = br
                }
            } else {
                if (i == 0) {
                    break
                }
                a = EDGE_STACK[--i]
            }
        }

        return ar
    }

    private fun inCircle(a: PointF, b: PointF, c: PointF, p: PointF): Boolean {

        // distance between a and p
        val dx = a.x - p.x
        val dy = a.y - p.y
        val ap = dx * dx + dy * dy

        // distance between b and p
        val ex = b.x - p.x
        val ey = b.y - p.y
        val bp = ex * ex + ey * ey

        // distance between c and p
        val fx = c.x - p.x
        val fy = c.y - p.y
        val cp = fx * fx + fy * fy

        return dx * (ey * cp - bp * fy) -
                dy * (ex * cp - bp * fx) +
                ap * (ex * fy - ey * fx) < 0
    }
}