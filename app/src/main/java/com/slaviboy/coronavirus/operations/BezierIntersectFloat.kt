package com.slaviboy.coronavirus.operations

// Copyright (C) 2020 Stanislav Georgiev
//  https://github.com/slaviboy
//
//	This program is free software: you can redistribute it and/or modify
//	it under the terms of the GNU Affero General Public License as
//	published by the Free Software Foundation, either version 3 of the
//	License, or (at your option) any later version.
//
//	This program is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU Affero General Public License for more details.
//
//	You should have received a copy of the GNU Affero General Public License
//	along with this program.  If not, see <http://www.gnu.org/licenses/>.


/**
 * Object for determining the points of intersection between a
 * cubic bezier curve and a line.
 */
object BezierIntersectFloat {

    var POLYNOMIAL_TOLERANCE = 1e-6
    var TOLERANCE = 1e-12

    private fun getPolynomialRoots(vararg C: Float): ArrayList<Float> {
        var degree = C.size - 1
        val n = degree
        val results = ArrayList<Float>()
        for (i in 0 until degree + 1) {
            if (Math.abs(C[i]) <= TOLERANCE) {
                degree--
            } else {
                break
            }
        }

        when (degree) {
            1 -> getLinearRoots(C[n], C[n - 1], results)
            2 -> getQuadraticRoots(C[n], C[n - 1], C[n - 2], results)
            3 -> getCubicRoots(C[n], C[n - 1], C[n - 2], C[n - 3], results)
        }

        return results
    }

    private fun getLinearRoots(
        C0: Float, C1: Float, results: ArrayList<Float> = ArrayList()
    ): ArrayList<Float> {

        if (C1 != 0.0f) {
            results.add(-C0 / C1)
        }
        return results
    }

    private fun getQuadraticRoots(
        C0: Float, C1: Float, C2: Float,
        results: ArrayList<Float> = ArrayList<Float>()
    ): ArrayList<Float> {

        val a = C2
        val b = C1 / a
        val c = C0 / a
        val d = b * b - 4.0f * c

        if (d > 0) {
            val e = Math.sqrt(d.toDouble()).toFloat()
            results.add(0.5f * (-b + e))
            results.add(0.5f * (-b - e))
        } else if (d == 0.0f) {
            results.add(0.5f * -b)
        }

        return results
    }

    private fun getCubicRoots(
        C0: Float, C1: Float, C2: Float, C3: Float,
        results: ArrayList<Float> = ArrayList<Float>()
    ): ArrayList<Float> {

        val c3 = C3
        val c2 = C2 / c3
        val c1 = C1 / c3
        val c0 = C0 / c3

        val a = (3.0f * c1 - c2 * c2) / 3.0f
        val b = (2.0f * c2 * c2 * c2 - 9.0f * c1 * c2 + 27.0f * c0) / 27.0f
        val offset = c2 / 3.0f
        var discrim = b * b / 4.0f + a * a * a / 27.0f
        val halfB = b / 2.0f
        var tmp = 0.0f
        var root = 0.0f

        if (Math.abs(discrim) <= POLYNOMIAL_TOLERANCE) {
            discrim = 0.0f
        }

        if (discrim > 0.0f) {
            val e = Math.sqrt(discrim.toDouble()).toFloat()

            tmp = -halfB + e
            if (tmp >= 0) {
                root = Math.pow(tmp.toDouble(), 1.0 / 3.0).toFloat()
            } else {
                root = -Math.pow(-tmp.toDouble(), 1.0 / 3.0).toFloat()
            }

            tmp = -halfB - e
            if (tmp >= 0) {
                root += Math.pow(tmp.toDouble(), 1.0 / 3.0).toFloat()
            } else {
                root -= Math.pow(-tmp.toDouble(), 1.0 / 3.0).toFloat()
            }

            results.add(root - offset)
        } else if (discrim < 0.0f) {
            val distance = Math.sqrt(-a / 3.0).toFloat()
            val angle = Math.atan2(Math.sqrt(-discrim.toDouble()), -halfB.toDouble()) / 3.0f
            val cos = Math.cos(angle).toFloat()
            val sin = Math.sin(angle).toFloat()
            val sqrt3 = Math.sqrt(3.0).toFloat()

            results.add(2.0f * distance * cos - offset)
            results.add(-distance * (cos + sqrt3 * sin) - offset)
            results.add(-distance * (cos - sqrt3 * sin) - offset)
        } else {
            if (halfB >= 0.0f) {
                tmp = -Math.pow(halfB.toDouble(), 1.0 / 3.0).toFloat()
            } else {
                tmp = Math.pow(-halfB.toDouble(), 1.0 / 3.0).toFloat()
            }

            results.add(2.0f * tmp - offset)
            // really should return next root twice, but we return only one
            results.add(-tmp - offset)
        }

        return results
    }

    fun quadBezierLine(
        p1x: Float, p1y: Float,
        p2x: Float, p2y: Float,
        p3x: Float, p3y: Float,
        a1x: Float, a1y: Float,
        a2x: Float, a2y: Float,
        result: ArrayList<Float>?
    ): Int {

        // temporary variables
        var ax = 0.0f
        var ay = 0.0f
        var bx = 0.0f
        var by = 0.0f

        // coefficients of quadratic
        var c2x = 0.0f
        var c2y = 0.0f
        var c1x = 0.0f
        var c1y = 0.0f
        var c0x = 0.0f
        var c0y = 0.0f

        // c coefficient for normal form of line
        var cl = 0.0f

        // normal for normal form of line
        var nx = 0.0f
        var ny = 0.0f

        // used to determine if point is on line segment
        val minx = Math.min(a1x, a2x)
        val miny = Math.min(a1y, a2y)
        val maxx = Math.max(a1x, a2x)
        val maxy = Math.max(a1y, a2y)

        ax = p2x * -2.0f
        ay = p2y * -2.0f
        c2x = p1x + ax + p3x
        c2y = p1y + ay + p3y

        ax = p1x * -2.0f
        ay = p1y * -2.0f
        bx = p2x * 2.0f
        by = p2y * 2.0f
        c1x = ax + bx
        c1y = ay + by

        // vec
        c0x = p1x
        c0y = p1y

        // Convert line to normal form: ax + by + c = 0
        // Find normal to line: negative inverse of original line's slope
        nx = a1y - a2y
        ny = a2x - a1x

        // Determine new c coefficient
        cl = a1x * a2y - a2x * a1y

        // Transform cubic coefficients to line's coordinate system
        // and find roots of cubic
        val roots = getPolynomialRoots(
            // dot products => x * x + y * y
            nx * c2x + ny * c2y,
            nx * c1x + ny * c1y,
            nx * c0x + ny * c0y + cl
        )

        // Any roots in closed interval [0,1] are intersections on Bezier, but
        // might not be on the line segment.
        // Find intersections and calculate point coordinates
        for (i in 0 until roots.size) {
            val t = roots[i]

            // We're within the Bezier curve
            if (0.0f <= t && t <= 1.0f) {

                // Find point on Bezier
                // lerp: x1 + (x2 - x1) * t
                val p4x = p1x + (p2x - p1x) * t
                val p4y = p1y + (p2y - p1y) * t

                val p5x = p2x + (p3x - p2x) * t
                val p5y = p2y + (p3y - p2y) * t

                // candidate
                val p6x = p4x + (p5x - p4x) * t
                val p6y = p4y + (p5y - p4y) * t

                // See if point is on line segment
                // Had to make special cases for vertical and horizontal lines due
                // to slight errors in calculation of p6
                if (a1x == a2x) {
                    if (miny <= p6y && p6y <= maxy) {
                        if (result != null) {
                            result.add(p6x)
                            result.add(p6y)
                        } else {
                            return 1
                        }
                    }
                } else if (a1y == a2y) {
                    if (minx <= p6x && p6x <= maxx) {
                        if (result != null) {
                            result.add(p6x)
                            result.add(p6y)
                        } else {
                            return 1
                        }
                    }

                    // gte: (x1 >= x2 && y1 >= y2)
                    // lte: (x1 <= x2 && y1 <= y2)
                } else if (p6x >= minx && p6y >= miny && p6x <= maxx && p6y <= maxy) {
                    if (result != null) {
                        result.add(p6x)
                        result.add(p6y)
                    } else {
                        return 1
                    }
                }
            }
        }
        return if (result != null) result.size / 2 else 0
    }

    fun quadBezierAABB(
        ax: Float, ay: Float,
        c1x: Float, c1y: Float,
        bx: Float, by: Float,
        xmin: Float, ymin: Float,
        xmax: Float, ymax: Float,
        result: ArrayList<Float>?
    ): Int {
        if (result != null) {
            // all intersections
            quadBezierLine(ax, ay, c1x, c1y, bx, by, xmin, ymin, xmax, ymin, result)
            quadBezierLine(ax, ay, c1x, c1y, bx, by, xmax, ymin, xmax, ymax, result)
            quadBezierLine(ax, ay, c1x, c1y, bx, by, xmin, ymax, xmax, ymax, result)
            quadBezierLine(ax, ay, c1x, c1y, bx, by, xmin, ymin, xmin, ymax, result)
            return result.size / 2
        } else {
            // any intersections
            // trivial cases
            if (xmin <= ax && xmax >= ax && ymin <= ay && ymax >= ay) {
                return 1
            }
            if (xmin <= bx && xmax >= bx && ymin <= by && ymax >= by) {
                return 1
            }
            if (quadBezierLine(ax, ay, c1x, c1y, bx, by, xmin, ymin, xmax, ymin, null) != 0) {
                return 1
            }
            if (quadBezierLine(ax, ay, c1x, c1y, bx, by, xmax, ymin, xmax, ymax, null) != 0) {
                return 1
            }
            if (quadBezierLine(ax, ay, c1x, c1y, bx, by, xmin, ymax, xmax, ymax, null) != 0) {
                return 1
            }
            if (quadBezierLine(ax, ay, c1x, c1y, bx, by, xmin, ymin, xmin, ymax, null) != 0) {
                return 1
            }
            return 0
        }
    }

    fun cubicBezierLine(
        p1x: Float, p1y: Float,
        p2x: Float, p2y: Float,
        p3x: Float, p3y: Float,
        p4x: Float, p4y: Float,
        a1x: Float, a1y: Float,
        a2x: Float, a2y: Float,
        result: ArrayList<Float>?
    ): Int {

        //temporary variables
        var ax = 0.0f
        var ay = 0.0f
        var bx = 0.0f
        var by = 0.0f
        var cx = 0.0f
        var cy = 0.0f
        var dx = 0.0f
        var dy = 0.0f

        // coefficients of cubic
        var c3x = 0.0f
        var c3y = 0.0f
        var c2x = 0.0f
        var c2y = 0.0f
        var c1x = 0.0f
        var c1y = 0.0f
        var c0x = 0.0f
        var c0y = 0.0f

        // c coefficient for normal form of line
        var cl = 0.0f

        // normal for normal form of line
        var nx = 0.0f
        var ny = 0.0f

        // used to determine if point is on line segment
        var minx = Math.min(a1x, a2x)
        var miny = Math.min(a1y, a2y)
        var maxx = Math.max(a1x, a2x)
        var maxy = Math.max(a1y, a2y)

        // Start with Bezier using Bernstein polynomials for weighting functions:
        //     (1-t^3)P1 + 3t(1-t)^2P2 + 3t^2(1-t)P3 + t^3P4
        //
        // Expand and collect terms to form linear combinations of original Bezier
        // controls.  This ends up with a vector cubic in t:
        //     (-P1+3P2-3P3+P4)t^3 + (3P1-6P2+3P3)t^2 + (-3P1+3P2)t + P1
        //             /\                  /\                /\       /\
        //             ||                  ||                ||       ||
        //             c3                  c2                c1       c0

        // Calculate the coefficients
        ax = p1x * -1.0f
        ay = p1y * -1.0f
        bx = p2x * 3.0f
        by = p2y * 3.0f
        cx = p3x * -3.0f
        cy = p3y * -3.0f
        dx = ax + bx + cx + p4x
        dy = ay + by + cy + p4y
        c3x = dx
        c3y = dy // vec

        ax = p1x * 3.0f
        ay = p1y * 3.0f
        bx = p2x * -6.0f
        by = p2y * -6.0f
        cx = p3x * 3.0f
        cy = p3y * 3.0f
        dx = ax + bx + cx
        dy = ay + by + cy
        c2x = dx
        c2y = dy // vec

        ax = p1x * -3.0f
        ay = p1y * -3.0f
        bx = p2x * 3.0f
        by = p2y * 3.0f
        cx = ax + bx
        cy = ay + by
        c1x = cx
        c1y = cy // vec

        c0x = p1x
        c0y = p1y

        // Convert line to normal form: ax + by + c = 0
        // Find normal to line: negative inverse of original line's slope
        nx = a1y - a2y
        ny = a2x - a1x

        // Determine new c coefficient
        cl = a1x * a2y - a2x * a1y

        // ?Rotate each cubic coefficient using line for new coordinate system?
        // Find roots of rotated cubic
        val roots = getPolynomialRoots(
            // dot products => x * x + y * y
            nx * c3x + ny * c3y,
            nx * c2x + ny * c2y,
            nx * c1x + ny * c1y,
            nx * c0x + ny * c0y + cl
        )

        // Any roots in closed interval [0,1] are intersections on Bezier, but
        // might not be on the line segment.
        // Find intersections and calculate point coordinates
        for (i in 0 until roots.size) {
            val t = roots[i]

            // We're within the Bezier curve
            if (0.0f <= t && t <= 1.0f) {

                // Find point on Bezier
                // lerp: x1 + (x2 - x1) * t
                val p5x = p1x + (p2x - p1x) * t
                val p5y = p1y + (p2y - p1y) * t // lerp(p1, p2, t);

                val p6x = p2x + (p3x - p2x) * t
                val p6y = p2y + (p3y - p2y) * t

                val p7x = p3x + (p4x - p3x) * t
                val p7y = p3y + (p4y - p3y) * t

                val p8x = p5x + (p6x - p5x) * t
                val p8y = p5y + (p6y - p5y) * t

                val p9x = p6x + (p7x - p6x) * t
                val p9y = p6y + (p7y - p6y) * t

                // candidate
                val p10x = p8x + (p9x - p8x) * t
                val p10y = p8y + (p9y - p8y) * t

                // See if point is on line segment
                if (a1x == a2x) {

                    // vertical
                    if (miny <= p10y && p10y <= maxy) {
                        if (result != null) {
                            result.add(p10x)
                            result.add(p10y)
                        } else {
                            return 1
                        }
                    }
                } else if (a1y == a2y) {

                    // horizontal
                    if (minx <= p10x && p10x <= maxx) {
                        if (result != null) {
                            result.add(p10x)
                            result.add(p10y)
                        } else {
                            return 1
                        }
                    }
                } else if (p10x >= minx && p10y >= miny && p10x <= maxx && p10y <= maxy) {
                    if (result != null) {
                        result.add(p10x)
                        result.add(p10y)
                    } else {
                        return 1
                    }
                }
            }
        }
        return if (result != null) result.size / 2 else 0
    }


    fun cubicBezierAABB(
        ax: Float, ay: Float,
        c1x: Float, c1y: Float,
        c2x: Float, c2y: Float,
        bx: Float, by: Float,
        xmin: Float, ymin: Float,
        xmax: Float, ymax: Float,
        result: ArrayList<Float>?
    ): Int {
        if (result != null) {

            // all intersections
            cubicBezierLine(ax, ay, c1x, c1y, c2x, c2y, bx, by, xmin, ymin, xmax, ymin, result)
            cubicBezierLine(ax, ay, c1x, c1y, c2x, c2y, bx, by, xmax, ymin, xmax, ymax, result)
            cubicBezierLine(ax, ay, c1x, c1y, c2x, c2y, bx, by, xmin, ymax, xmax, ymax, result)
            cubicBezierLine(ax, ay, c1x, c1y, c2x, c2y, bx, by, xmin, ymin, xmin, ymax, result)
            return result.size / 2
        } else {
            // any intersections
            // trivial cases
            if (xmin <= ax && xmax >= ax && ymin <= ay && ymax >= ay) {
                return 1
            }
            if (xmin <= bx && xmax >= bx && ymin <= by && ymax >= by) {
                return 1
            }
            if (cubicBezierLine(
                    ax,
                    ay,
                    c1x,
                    c1y,
                    c2x,
                    c2y,
                    bx,
                    by,
                    xmin,
                    ymin,
                    xmax,
                    ymin,
                    null
                ) != 0
            ) {
                return 1
            }
            if (cubicBezierLine(
                    ax,
                    ay,
                    c1x,
                    c1y,
                    c2x,
                    c2y,
                    bx,
                    by,
                    xmax,
                    ymin,
                    xmax,
                    ymax,
                    null
                ) != 0
            ) {
                return 1
            }
            if (cubicBezierLine(
                    ax,
                    ay,
                    c1x,
                    c1y,
                    c2x,
                    c2y,
                    bx,
                    by,
                    xmin,
                    ymax,
                    xmax,
                    ymax,
                    null
                ) != 0
            ) {
                return 1
            }
            if (cubicBezierLine(
                    ax,
                    ay,
                    c1x,
                    c1y,
                    c2x,
                    c2y,
                    bx,
                    by,
                    xmin,
                    ymin,
                    xmin,
                    ymax,
                    null
                ) != 0
            ) {
                return 1
            }
            return 0;
        }
    }

}