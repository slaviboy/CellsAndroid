package com.slaviboy.voronoifloat

import android.graphics.PointF
import android.graphics.RectF
import java.lang.Exception
import kotlin.collections.ArrayList

object Polygon {

    /**
     * Calculate the area of the polygon
     */
    fun area(points: ArrayList<PointF>): Float {

        var n = points.size
        var i = 0
        var a: PointF
        var b = points[n - 1]
        var area = 0.0

        while (i < n) {
            a = b
            b = points[i]
            area += a.y * b.x - a.x * b.y
            i++
        }
        return Math.abs(area / 2.0).toFloat()
    }

    /**
     * Returns the centroid point of the polygon
     */
    fun centroid(points: ArrayList<PointF>): PointF {

        val n = points.size
        var i = 0
        var x = 0.0f
        var y = 0.0f
        var a: PointF
        var b = points[n - 1]
        var c: Float
        var k = 0.0f

        while (i < n) {
            a = b
            b = points[i]
            c = a.x * b.y - b.x * a.y
            k += c
            x += (a.x + b.x) * c
            y += (a.y + b.y) * c
            i++
        }
        k *= 3.0f
        return PointF(x / k, y / k)
    }

    fun contains(points: ArrayList<PointF>, point: PointF): Boolean {
        var n = points.size
        var p = points[n - 1]
        var x = point.x
        var y = point.y
        var x0 = p.x
        var y0 = p.y
        var x1: Float
        var y1: Float
        var inside = false

        for (i in 0 until n) {
            var p = points[i]
            x1 = p.x
            y1 = p.y
            if (((y1 > y) != (y0 > y)) &&
                (x < (x0 - x1) * (y - y1) / (y0 - y1) + x1)
            ) {
                inside = !inside
            }
            x0 = x1
            y0 = y1
        }
        return inside
    }

    fun length(points: ArrayList<PointF>): Float {
        var i = -1
        var n = points.size
        var b = points[n - 1]
        var xa: Float
        var ya: Float
        var xb = b.x
        var yb = b.y
        var perimeter = 0.0f

        while (++i < n) {
            xa = xb
            ya = yb
            b = points[i]
            xb = b.x
            yb = b.y
            xa -= xb
            ya -= yb
            perimeter += Math.sqrt(0.0 + xa * xa + ya * ya).toFloat()
        }
        return perimeter
    }

    // Returns the 2D cross product of AB and AC vectors, i.e., the z-component of
    // the 3D cross product in a quadrant I Cartesian coordinate system (+x is
    // right, +y is up). Returns a positive value if ABC is counter-clockwise,
    // negative if clockwise, and zero if the points are collinear.
    fun cross(a: SortedPointF, b: SortedPointF, c: SortedPointF): Float {
        return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
    }

    // Computes the upper convex hull per the monotone chain algorithm.
    // Assumes points.length >= 3, is sorted by x, unique in y.
    // Returns an array of indices into points in left-to-right order.
    fun computeUpperHullIndexes(points: Array<SortedPointF>): ArrayList<Int> {

        var n = points.size
        val indexes = ArrayList<Int>()
        indexes.add(0)
        indexes.add(1)
        var size = 2

        for (i in 2 until n) {
            while (size > 1 &&
                cross(
                    points[indexes[size - 2]],
                    points[indexes[size - 1]],
                    points[i]
                ) <= 0
            ) {
                size--
            }

            indexes.add(size, i)
            size++
        }

        return indexes.slice(0 until size) as ArrayList<Int> // remove popped points
    }

    fun hull(points: ArrayList<PointF>): ArrayList<PointF> {

        val n = points.size
        if (n < 3) throw Exception("Polygon must have at least 3 points")

        val sortedPoints: Array<SortedPointF> = Array(n) {
            SortedPointF(points[it].x, points[it].y, it)
        }
        sortedPoints.sortedWith(compareBy({ it.x }, { it.y }))

        val flippedPoints: Array<SortedPointF> = Array(n) {
            SortedPointF(sortedPoints[it].x, -sortedPoints[it].y, 0)
        }

        val upperIndexes = computeUpperHullIndexes(sortedPoints)
        val lowerIndexes = computeUpperHullIndexes(flippedPoints)

        // Construct the hull polygon, removing possible duplicate endpoints.
        val skipLeft = if (lowerIndexes[0] == upperIndexes[0]) 1 else 0
        val skipRight =
            if (lowerIndexes[lowerIndexes.size - 1] == upperIndexes[upperIndexes.size - 1]) 1 else 0
        val hull = ArrayList<PointF>()

        // Add upper hull in right-to-l order.
        // Then add lower hull in left-to-right order.
        for (i in upperIndexes.indices.reversed()) {
            hull.add(points[sortedPoints[upperIndexes[i]]!!.index])
        }
        for (i in skipLeft until lowerIndexes.size - skipRight) {
            hull.add(points[sortedPoints[lowerIndexes[i]]!!.index])
        }

        return hull
    }


    /**
     * Get the bound rectangle of a polygon, and return it top, left, right and bottom position
     * @param points polygon points
     */
    fun bound(points: ArrayList<PointF>): RectF {
        var left = Float.MAX_VALUE
        var rights = Float.MIN_VALUE
        var top = Float.MAX_VALUE
        var bottom = Float.MIN_VALUE

        points.forEach {
            if (it.x < left) {
                left = it.x
            }

            if (it.x > rights) {
                rights = it.x
            }

            if (it.y < top) {
                top = it.y
            }

            if (it.y > bottom) {
                bottom = it.y
            }
        }

        return RectF(left, top, rights, bottom)
    }

}

class SortedPointF(var x: Float = 0f, var y: Float = 0f, var index: Int = 0)// : PointF(x, y)