package com.slaviboy.coronavirus.operations

import android.graphics.*
import com.slaviboy.voronoifloat.Polygon
import com.slaviboy.voronoifloat.Voronoi
import java.lang.Exception

class LineF(var a: PointF = PointF(), var b: PointF = PointF(), var index: Int = 0)

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
 * Class that generates mesh bitmap, for particular bitmap to fit inside a
 * polygon path with rounded corners. The rounded path is made out of cubic
 * bezier curves.
 *
 * @param bitmap bitmaps, that will be the texture for the cells
 * @param pointsPerWidth how many mesh points for the bitmap width(rows)
 * @param pointsPerHeight how many mesh points for the bitmap height(column)
 * @param voronoi the voronoi object
 * @param scale scale the paths, as the pivot point is the center of the path
 */
class Cells(
    var bitmap: Bitmap,
    var pointsPerWidth: Int,
    var pointsPerHeight: Int,
    var voronoi: Voronoi,
    var voronoiStrokeColor: Int = Color.BLACK,
    var voronoiFillColor: Int = Color.TRANSPARENT,
    var curveVoronoiStrokeColor: Int = Color.RED,
    var curveVoronoiFillColor: Int = Color.TRANSPARENT,
    var scale: PointF = PointF(0.9f, 0.9f),
    var cubicCurvePointStrokeColor: Int = Color.TRANSPARENT,
    var cubicCurvePointFillColor: Int = Color.BLUE,
    var cubicCurvePointRadius: Float = 8.0f,
    var centerPointStrokeColor: Int = Color.TRANSPARENT,
    var centerPointFillColor: Int = Color.RED,
    var centerPointRadius: Float = 8.0f
) {

    private val totalPoints = (pointsPerWidth + 1) * (pointsPerHeight + 1)
    private val originalMeshPoints: Array<PointF> = Array(totalPoints) { PointF() }
    private var lines: ArrayList<ArrayList<LineF>>
    private var indicesPerLayer: ArrayList<ArrayList<Int>>
    private var indicesLayerTopRightBottomLeft: ArrayList<ArrayList<Int>>

    init {
        generateOriginalMeshPoints()

        // get all lines from the center to each grid point for the bitmap
        lines = getCenterToPointLines(
            PointF(bitmap.width / 2.0f, bitmap.height / 2.0f),
            originalMeshPoints, pointsPerWidth + 1, pointsPerHeight + 1
        )

        indicesPerLayer = getIndicesPerLayer(pointsPerWidth + 1, pointsPerHeight + 1)
        indicesLayerTopRightBottomLeft =
            getIndicesLayerTopRightBottomLeft(pointsPerWidth + 1, pointsPerHeight + 1)
    }

    /**
     * Generate the original mesh points, for the number of
     * pointsPerWidth and pointsPerHeight.
     */
    private fun generateOriginalMeshPoints() {

        val bitmapWidth: Float = bitmap.width.toFloat()
        val bitmapHeight: Float = bitmap.height.toFloat()

        var index = 0
        for (y in 0..pointsPerHeight) {
            val fy = bitmapHeight * y / pointsPerHeight
            for (x in 0..pointsPerWidth) {
                val fx = bitmapWidth * x / pointsPerWidth
                originalMeshPoints[index] = PointF(fx, fy)
                index++
            }
        }
    }

    /**
     * Returns double array with indices from each layer
     */
    fun getIndicesPerLayer(rows: Int, columns: Int): ArrayList<ArrayList<Int>> {

        val hS = 0
        val hE = columns

        val wS = 0
        val wE = rows

        if (hE - hS < 0 || wE - wS < 0) {
            throw Exception("Not enough elements !!!")
        }

        val pointsIndex = ArrayList<ArrayList<Int>>()
        val firstLayer = ArrayList<Int>()

        // top
        for (i in wS until wE) {
            val pointIndex = rows * hS + i
            firstLayer.add(pointIndex)
        }

        // right
        for (i in hS + 1 until hE) {
            val pointIndex = rows * i + (wE - 1)
            firstLayer.add(pointIndex)
        }

        // bottom
        for (i in wE - 2 downTo wS) {
            val pointIndex = rows * (hE - 1) + i
            firstLayer.add(pointIndex)
        }

        // left
        for (i in hE - 2 downTo hS + 1) {
            val pointIndex = rows * i + wS
            firstLayer.add(pointIndex)
        }

        pointsIndex.add(firstLayer)

        val min = Math.min(rows, columns)

        var previousLayer = firstLayer
        val addOdd = 2

        for (i in 1 until min / 2) {

            val addEven = min - 3 - (2 * (i - 1))
            val layerRemove = ArrayList<Int>()
            var sum = 1 - 2
            for (j in 0 until 8) {
                val h = if (j % 2 == 0) addOdd else addEven
                sum += h
                layerRemove.add(previousLayer[sum])
            }

            val newLayer = (previousLayer.clone() as ArrayList<Int>)
            newLayer.removeAll(layerRemove)

            pointsIndex.add(newLayer)

            previousLayer = newLayer
        }

        return pointsIndex
    }

    /**
     * Returns double array with indices for each layer, that correspond to indices
     * from the fist layer with remove indices after the path coordinates are scaled down.
     * @param rows number of row points
     * @param columns number of columns points
     */
    fun getIndicesLayerTopRightBottomLeft(rows: Int, columns: Int): ArrayList<ArrayList<Int>> {

        val indices = ArrayList<ArrayList<Int>>()
        var index = 0
        while (true) {

            val hS = index
            val hE = columns - 1 - index

            val wS = index
            val wE = rows - 1 - index

            if (hE - hS < 0 || wE - wS < 0) {
                break
            }

            val indicesLayer = ArrayList<Int>()

            // top
            for (i in wS until wE + 1) {
                val pointIndex = rows * hS + i
                indicesLayer.add(pointIndex)
            }

            // right
            for (i in hS + 1 until hE + 1) {
                val pointIndex = rows * i + wE
                indicesLayer.add(pointIndex)
            }

            // bottom
            for (i in wE - 1 downTo wS) {
                val pointIndex = rows * hE + i
                indicesLayer.add(pointIndex)
            }

            // left
            for (i in hE - 1 downTo hS + 1) {
                val pointIndex = rows * i + wS
                indicesLayer.add(pointIndex)
            }

            indices.add(indicesLayer)

            index++
        }
        return indices
    }

    /**
     * Matrix map method for array list with floating points, it applies the transformation
     * from the matrix directly to the same array with points.
     *
     * @param points array list with points that will be transformed
     */
    fun Matrix.mapPoints(points: ArrayList<PointF>) {

        val inPoint = FloatArray(points.size * 2)
        for (i in 0 until points.size) {
            inPoint[i * 2] = points[i].x
            inPoint[i * 2 + 1] = points[i].y
        }
        this.mapPoints(inPoint)

        for (i in 0 until points.size) {
            points[i] = PointF(
                inPoint[i * 2],
                inPoint[i * 2 + 1]
            )
        }
    }

    /**
     * Matrix map method applying transformation from the matrix to the transformedPoints
     * array list, for the points from the array list -points.
     *
     * @param points array list with points that will be transformed
     * @param transformedPoints array list with the transformed points
     */
    fun Matrix.mapPoints(points: ArrayList<PointF>, transformedPoints: ArrayList<PointF>) {

        val inPoints = FloatArray(points.size * 2)
        for (i in 0 until points.size) {
            inPoints[i * 2] = points[i].x
            inPoints[i * 2 + 1] = points[i].y
        }
        this.mapPoints(inPoints)

        for (i in 0 until points.size) {
            val point = PointF(
                inPoints[i * 2],
                inPoints[i * 2 + 1]
            )
            transformedPoints.add(point)
        }
    }

    fun generateCubicCurve(cubicCurvePoints: ArrayList<PointF>): ArrayList<ArrayList<PointF>> {

        // filter
        val filterPoints = ArrayList<PointF>()
        for (i in 0 until cubicCurvePoints.size) {
            val j = (i + 1) % cubicCurvePoints.size
            val d = distance(cubicCurvePoints[i], cubicCurvePoints[j]);
            if (d < 10) {
                cubicCurvePoints[j] = PointF(
                    0.5f * cubicCurvePoints[i].x + 0.5f * cubicCurvePoints[j].x,
                    0.5f * cubicCurvePoints[i].y + 0.5f * cubicCurvePoints[j].y
                )
            } else {
                filterPoints.add(PointF(cubicCurvePoints[i].x, cubicCurvePoints[i].y))
            }
        }

        val middlePoints = ArrayList<PointF>()
        for (i in 0 until filterPoints.size) {
            val j = (i + 1) % filterPoints.size
            middlePoints.add(
                PointF(
                    0.5f * filterPoints[i].x + 0.5f * filterPoints[j].x,
                    0.5f * filterPoints[i].y + 0.5f * filterPoints[j].y
                )
            )
        }

        // add coordinates for each cubic bezier curve
        val cubicCurves = ArrayList<ArrayList<PointF>>()
        var previousPoint = middlePoints[0]
        for (i in 0 until middlePoints.size) {
            val j = (i + 1) % middlePoints.size
            cubicCurves.add(
                arrayListOf(
                    previousPoint,
                    filterPoints[j],
                    filterPoints[j],
                    middlePoints[j]
                )
            )
            previousPoint = middlePoints[j]
        }

        return cubicCurves
    }

    /**
     * Distance between two points
     * @param a first point
     * @param b second point
     */
    fun distance(a: PointF, b: PointF): Float {
        return Math.sqrt(0.0 + (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y))
            .toFloat()
    }

    /**
     * Get the new mesh coordinates for transforming a bitmap via the canvas
     * method drawBitmapMesh(). Those coordinate are transformed to fit inside the
     * curved path, formed by the polygon points from the voronoi cell.
     *
     * @param points polygon points from the voronoi cell
     * !!!slower algorithm
     */
    fun getMeshCoordinates(points: ArrayList<PointF>): FloatArray {
        return drawMeshCoordinates(points = points)
    }

    /**
     * Get the new mesh coordinates for transforming a bitmap via the canvas
     * method drawBitmapMesh(). Those coordinate are transformed to fit inside the
     * curved path, formed by the polygon points from the voronoi cell.
     *
     * @param points polygon points from the voronoi cell
     */
    fun getMeshCoordinatesFaster(points: ArrayList<PointF>): FloatArray {
        return drawMeshCoordinatesFaster(points = points)
    }

    /**
     * Faster algorithm for generating the mesh coordinates, and
     * drawing the generated bitmap on the canvas.
     */
    fun drawMeshCoordinatesFaster(canvas: Canvas? = null, points: ArrayList<PointF>): FloatArray {

        val cubicCurve = generateCubicCurve(points)

        var bound = Polygon.bound(points)
        val boundWidth = bound.width()
        val boundHeight = bound.height()
        val widthIsBigger = boundWidth > boundHeight
        val max = Math.max(boundWidth, boundHeight)

        var totalScales = 0
        var scale = 0.0f
        if (widthIsBigger) {
            scale = bitmap.width / max
            totalScales = (pointsPerWidth) / 2
        } else {
            scale = bitmap.height / max
            totalScales = (pointsPerHeight) / 2
        }

        // transform matrix to fit to bitmap size and translate to its center
        val matrix = Matrix().apply {
            postTranslate(-bound.left, -bound.top)
            postScale(scale, scale, bound.centerX() - bound.left, bound.centerY() - bound.top)
            postTranslate(
                (bitmap.width - bound.width()) / 2f,
                (bitmap.height - bound.height()) / 2f
            )
        }

        val matrixInvert = Matrix()
        matrix.invert(matrixInvert)

        // transform all points from cubic curve
        cubicCurve.forEach {
            matrix.mapPoints(it)
        }

        // get the new points to calculate the bound
        val pointsNew = ArrayList<PointF>()
        matrix.mapPoints(points, pointsNew)

        // calculate the new bound
        bound = Polygon.bound(pointsNew)

        val meshCoordinates = FloatArray(totalPoints * 2)
        val meshCoordinatesFitCurve = FloatArray(totalPoints * 2)
        for (i in originalMeshPoints.indices) {
            meshCoordinates[i * 2] = originalMeshPoints[i].x
            meshCoordinates[i * 2 + 1] = originalMeshPoints[i].y

            meshCoordinatesFitCurve[i * 2] = originalMeshPoints[i].x
            meshCoordinatesFitCurve[i * 2 + 1] = originalMeshPoints[i].y
        }

        cubicCurve.forEach {

            for (i in 0 until indicesPerLayer[0].size) {

                val index = indicesPerLayer[0][i]
                val result = ArrayList<Float>()
                BezierIntersectFloat.cubicBezierLine(
                    it[0].x, it[0].y,
                    it[1].x, it[1].y,
                    it[2].x, it[2].y,
                    it[3].x, it[3].y,
                    bitmap.width / 2f, bitmap.height / 2f,
                    originalMeshPoints[index].x + 0.1f, originalMeshPoints[index].y + 0.1f,
                    result
                )

                if (result.size >= 2) {
                    meshCoordinatesFitCurve[index * 2] = result[0]
                    meshCoordinatesFitCurve[index * 2 + 1] = result[1]
                }
            }
        }


        for (k in 0 until totalScales  ) {

            val newScaleX = 1 - k / totalScales.toFloat()
            val newScaleY = 1 - k / totalScales.toFloat()

            matrix.setScale(newScaleX, newScaleY, bound.centerX(), bound.centerY())

            val newCoordinate = FloatArray(meshCoordinatesFitCurve.size)
            matrix.mapPoints(newCoordinate, meshCoordinatesFitCurve)

            for (i in indicesPerLayer[k].indices) {
                val index = indicesPerLayer[k][i]
                val index2 = indicesLayerTopRightBottomLeft[k][i]
                meshCoordinates[index2 * 2] = newCoordinate[index * 2]
                meshCoordinates[index2 * 2 + 1] = newCoordinate[index * 2 + 1]
            }
        }

        matrixInvert.mapPoints(meshCoordinates)
        canvas?.drawBitmapMesh(
            bitmap,
            pointsPerWidth,
            pointsPerHeight,
            meshCoordinates,
            0, null,
            0, null
        )

        return meshCoordinates
    }

    fun drawMeshCoordinates(
        canvas: Canvas? = null,
        points: ArrayList<PointF>
    ): FloatArray {

        val cubicCurve = generateCubicCurve(points)

        var bound = Polygon.bound(points)
        val boundWidth = bound.width()
        val boundHeight = bound.height()
        val widthIsBigger = boundWidth > boundHeight
        val max = Math.max(boundWidth, boundHeight)

        var totalScales = 0
        var scale = 0.0f
        if (widthIsBigger) {
            scale = bitmap.width / max
            totalScales = (pointsPerWidth) / 2
        } else {
            scale = bitmap.height / max
            totalScales = (pointsPerHeight) / 2
        }

        // transform matrix to fit to bitmap size and translate to its center
        val matrix = Matrix().apply {
            postTranslate(-bound.left, -bound.top)
            postScale(scale, scale, bound.centerX() - bound.left, bound.centerY() - bound.top)
            postTranslate(
                (bitmap.width - bound.width()) / 2f,
                (bitmap.height - bound.height()) / 2f
            )
        }

        val matrixInvert = Matrix()
        matrix.invert(matrixInvert)

        // transform all points from cubic curve
        cubicCurve.forEach {
            matrix.mapPoints(it)
        }

        // get the new points to calculate the bound
        val pointsNew = ArrayList<PointF>()
        matrix.mapPoints(points, pointsNew)

        // calculate the new bound
        bound = Polygon.bound(pointsNew)

        // get the center coordinate of the bound box
        val centerX = bound.centerX()
        val centerY = bound.centerY()

        val meshCoordinates = FloatArray(totalPoints * 2)
        for (i in originalMeshPoints.indices) {
            meshCoordinates[i * 2] = originalMeshPoints[i].x
            meshCoordinates[i * 2 + 1] = originalMeshPoints[i].y
        }

        for (k in 0 until totalScales) {

            val newScaleX = 1 - k / totalScales.toFloat()
            val newScaleY = 1 - k / totalScales.toFloat()

            matrix.setScale(newScaleX, newScaleY, centerX, centerY)

            val cubicCurveNew = ArrayList<ArrayList<PointF>>()
            for (i in cubicCurve.indices) {
                val newPoints = ArrayList<PointF>()
                matrix.mapPoints(cubicCurve[i], newPoints)
                cubicCurveNew.add(newPoints)
            }

            cubicCurveNew.forEach {
                if (k < lines.size) {
                    for (i in 0 until lines[k].size) {
                        val line = lines[k][i]

                        val result = ArrayList<Float>()
                        BezierIntersectFloat.cubicBezierLine(
                            it[0].x, it[0].y,
                            it[1].x, it[1].y,
                            it[2].x, it[2].y,
                            it[3].x, it[3].y,
                            line.a.x, line.a.y,
                            line.b.x + 0.1f, line.b.y + 0.1f,
                            result
                        )

                        if (result.size >= 2) {
                            meshCoordinates[line.index * 2] = result[0]
                            meshCoordinates[line.index * 2 + 1] = result[1]
                        }
                    }
                }
            }
        }

        matrixInvert.mapPoints(meshCoordinates)
        canvas?.drawBitmapMesh(
            bitmap,
            pointsPerWidth,
            pointsPerHeight,
            meshCoordinates,
            0, null,
            0, null
        )

        return meshCoordinates
    }

    fun drawMeshBitmap(canvas: Canvas, paint: Paint, meshCoordinates: FloatArray) {

        canvas.drawBitmapMesh(
            bitmap,
            pointsPerWidth,
            pointsPerHeight,
            meshCoordinates,
            0, null,
            0, paint
        )
    }


    /**
     * Get all lines from the center to the original mesh points
     * @param center center point of the bitmap
     * @param points array with all points from
     * @param rows number of points per row
     * @param columns number of points per column
     */
    fun getCenterToPointLines(
        center: PointF, points: Array<PointF>, rows: Int, columns: Int
    ): ArrayList<ArrayList<LineF>> {

        val sortedLines = ArrayList<ArrayList<LineF>>()
        var index = 0
        while (true) {

            val hS = index
            val hE = columns - 1 - index

            val wS = index
            val wE = rows - 1 - index

            if (hE - hS < 0 || wE - wS < 0) {
                break
            }

            val lines = ArrayList<LineF>()

            //FIX REPLACE PointF(point.x, point.y) with point !!
            // top
            for (i in wS until wE + 1) {
                val pointIndex = rows * hS + i
                val point = points[pointIndex]
                val line = LineF(center, PointF(point.x, point.y), pointIndex)
                lines.add(line)
            }

            // bottom
            for (i in wS until wE + 1) {
                val pointIndex = rows * hE + i
                val point = points[pointIndex]
                val line = LineF(center, PointF(point.x, point.y), pointIndex)
                lines.add(line)
            }

            // left
            for (i in hS until hE + 1) {
                val pointIndex = rows * i + wS
                val point = points[pointIndex]
                val line = LineF(center, PointF(point.x, point.y), pointIndex)
                lines.add(line)
            }

            // right
            for (i in hS until hE + 1) {
                val pointIndex = rows * i + wE
                val point = points[pointIndex]
                val line = LineF(center, PointF(point.x, point.y), pointIndex)
                lines.add(line)
            }

            sortedLines.add(lines)

            index++
        }
        return sortedLines
    }


    /**
     * Draw all voronoi cells
     */
    fun drawVoronoi(canvas: Canvas, paint: Paint, path: Path) {
        path.reset()
        voronoi.render(path)
        drawPath(voronoiFillColor, voronoiStrokeColor, canvas, paint, path)
    }

    /**
     * Draw voronoi with smooth curves, that are scaled, depending on the scale
     * factor from the properties object.
     */
    fun drawCurveVoronoi(canvas: Canvas, paint: Paint, path: Path) {

        path.reset()
        for (k in voronoi.delaunay.points.indices) {
            val cubicCurvePoints = voronoi.cell(k)
            val center = Polygon.centroid(cubicCurvePoints)

            val matrix = Matrix().apply {
                setScale(
                    scale.x,
                    scale.y,
                    center.x,
                    center.y
                )
            }

            // filter
            val filterPoints = ArrayList<PointF>()
            for (i in 0 until cubicCurvePoints.size) {
                val j = (i + 1) % cubicCurvePoints.size
                val d = distance(cubicCurvePoints[i], cubicCurvePoints[j]);
                if (d < 10) {
                    cubicCurvePoints[j] = PointF(
                        0.5f * cubicCurvePoints[i].x + 0.5f * cubicCurvePoints[j].x,
                        0.5f * cubicCurvePoints[i].y + 0.5f * cubicCurvePoints[j].y
                    )
                } else {
                    filterPoints.add(PointF(cubicCurvePoints[i].x, cubicCurvePoints[i].y))
                }
            }


            // scale points using matrix
            val inPoint = FloatArray(filterPoints.size * 2)
            for (i in 0 until filterPoints.size) {
                inPoint[i * 2] = filterPoints[i].x
                inPoint[i * 2 + 1] = filterPoints[i].y
            }
            matrix.mapPoints(inPoint)
            for (i in 0 until filterPoints.size) {
                filterPoints[i] = PointF(
                    inPoint[i * 2],
                    inPoint[i * 2 + 1]
                )
            }


            val middlePoints = ArrayList<PointF>()
            for (i in 0 until filterPoints.size) {
                val j = (i + 1) % filterPoints.size
                middlePoints.add(
                    PointF(
                        0.5f * filterPoints[i].x + 0.5f * filterPoints[j].x,
                        0.5f * filterPoints[i].y + 0.5f * filterPoints[j].y
                    )
                )
            }

            path.reset()
            path.moveTo(middlePoints[0].x, middlePoints[0].y)
            for (i in 0 until middlePoints.size) {
                val j = (i + 1) % middlePoints.size
                path.cubicTo(
                    filterPoints[j].x,
                    filterPoints[j].y,
                    filterPoints[j].x,
                    filterPoints[j].y,
                    middlePoints[j].x,
                    middlePoints[j].y
                )
            }

            drawPath(
                curveVoronoiFillColor,
                curveVoronoiStrokeColor,
                canvas,
                paint,
                path
            )
        }
    }

    /**
     * Draw all voronoi centers to all voronoi cells.
     */
    fun drawVoronoiCenters(canvas: Canvas, paint: Paint, path: Path) {
        path.reset()
        for (i in voronoi.delaunay.points.indices) {
            val cubicCurvePoints = voronoi.cell(i)
            val center = Polygon.centroid(cubicCurvePoints)
            path.addCircle(
                center.x,
                center.y,
                centerPointRadius,
                Path.Direction.CW
            )
        }
        drawPath(
            centerPointFillColor,
            centerPointStrokeColor,
            canvas,
            paint,
            path
        )
    }

    /**
     * Draw the cubicCurve voronoi points that were passed as arguments to the
     * Delaunay class, and are used to create the voronoi diagram.
     */
    fun drawVoronoiPoints(canvas: Canvas, paint: Paint, path: Path) {
        path.reset()
        for (i in voronoi.delaunay.points.indices) {
            val point = voronoi.delaunay.points[i]
            path.addCircle(
                point.x,
                point.y,
                cubicCurvePointRadius,
                Path.Direction.CW
            )
        }
        drawPath(
            cubicCurvePointFillColor,
            cubicCurvePointStrokeColor,
            canvas,
            paint,
            path
        )
    }

    /**
     * Draw fill and stroke path
     */
    fun drawPath(fillColor: Int, strokeColor: Int, canvas: Canvas, paint: Paint, path: Path) {

        if (fillColor != Color.TRANSPARENT) {
            paint.color = fillColor
            paint.style = Paint.Style.FILL
            canvas.drawPath(path, paint)
        }

        if (strokeColor != Color.TRANSPARENT) {
            paint.color = strokeColor
            paint.style = Paint.Style.STROKE
            canvas.drawPath(path, paint)
        }
    }


    fun fingerOverCell(x: Float, y: Float): PointF? {
        val points = voronoi.delaunay.points
        for (i in points.indices) {
            val point = points[i]
            val dx = Math.abs(x - point.x)
            val dy = Math.abs(y - point.y)
            if (dx < cubicCurvePointRadius * 4 && dy < cubicCurvePointRadius * 4) {
                return point
            }
        }
        return null
    }

    fun fingerOverCellIndex(x: Float, y: Float): Int {
        val points = voronoi.delaunay.points
        for (i in points.indices) {
            val point = points[i]
            val dx = Math.abs(x - point.x)
            val dy = Math.abs(y - point.y)
            if (dx < cubicCurvePointRadius * 4 && dy < cubicCurvePointRadius * 4) {
                return i
            }
        }
        return -1
    }
}