package com.slaviboy.coronavirus.drawing

import android.graphics.*
import android.view.MotionEvent
import android.view.View
import com.slaviboy.coronavirus.operations.Cells
import com.slaviboy.voronoifloat.Polygon

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
 * Simple drawing thread, used for drawing the particles on canvas specified by the
 * surface holder
 */
class DrawingThread(
    var surfaceHolder: ISurfaceHolder,
    var cells: Cells
) : Thread(), View.OnTouchListener {

    var isRunning = false
    private var previousTime: Long
    private val fps = 30

    var paint: Paint = Paint().apply {
        isAntiAlias = true
        color = Color.BLACK
        style = Paint.Style.STROKE
    }

    var lastPoints: ArrayList<PointF>
    var lastMeshCoordinates: ArrayList<FloatArray>
    var selectedPointIndex: Int = -1

    init {
        previousTime = System.currentTimeMillis()

        lastPoints = ArrayList()
        for (i in cells.voronoi.delaunay.points.indices) {
            lastPoints.add(PointF())
        }

        lastMeshCoordinates = ArrayList()
        for (i in cells.voronoi.delaunay.points.indices) {
            lastMeshCoordinates.add(FloatArray((cells.pointsPerWidth + 1) * (cells.pointsPerHeight + 1)))
        }
    }

    /**
     * Runnable method for the thread
     */
    override fun run() {
        var canvas: Canvas?
        while (isRunning) {

            val currentTimeMillis = System.currentTimeMillis()
            val elapsedTimeMs = currentTimeMillis - previousTime
            val sleepTimeMs = (1000f / fps - elapsedTimeMs).toLong()
            canvas = null

            try {
                canvas = surfaceHolder.lockCanvas()

                // sleep the thread if frame rate is bigger than the fps
                if (canvas == null) {
                    sleep(1)
                    continue
                } else if (sleepTimeMs > 0) {
                    sleep(sleepTimeMs)
                }

                // draw the particles
                synchronized(surfaceHolder) {

                    canvas.drawColor(Color.BLACK)
                    cells.voronoi.update()

                    for (i in cells.voronoi.delaunay.points.indices) {
                        val points = cells.voronoi.cell(i)
                        val center = Polygon.centroid(points)
                        if (lastPoints[i] != center) {

                            lastMeshCoordinates[i] = cells.getMeshCoordinatesFaster(points)
                            lastPoints[i] = center
                        }
                        cells.drawMeshBitmap(canvas, paint, lastMeshCoordinates[i])


                        canvas.drawCircle(center.x,center.y, cells.centerPointRadius, paint)
                    }
                }

            } catch (e: Exception) {
                e.printStackTrace()
            } finally {
                if (canvas != null) {
                    surfaceHolder.unlockCanvasAndPost(canvas)
                    previousTime = System.currentTimeMillis()
                }
            }
        }
    }

    override fun onTouch(view: View, motionEvent: MotionEvent): Boolean {

        when (motionEvent.action) {

            MotionEvent.ACTION_DOWN -> {
                selectedPointIndex = cells.fingerOverCellIndex(motionEvent.x, motionEvent.y)
            }

            MotionEvent.ACTION_UP -> {
                selectedPointIndex = -1
            }

            MotionEvent.ACTION_MOVE -> {
                if (selectedPointIndex != -1) {

                    // move the center points, for the cell over the finger
                    cells.voronoi.delaunay.points[selectedPointIndex] =
                        PointF(motionEvent.x, motionEvent.y)
                }
            }
        }

        return true
    }


}

