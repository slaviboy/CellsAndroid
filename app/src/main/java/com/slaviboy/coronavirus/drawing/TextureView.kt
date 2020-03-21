package com.slaviboy.coronavirus.drawing

import android.content.Context
import android.graphics.*
import android.util.AttributeSet
import android.view.TextureView
import android.view.TextureView.SurfaceTextureListener
import com.slaviboy.coronavirus.operations.Cells
import com.slaviboy.coronavirus.R
import com.slaviboy.voronoifloat.Delaunay
import com.slaviboy.voronoifloat.Voronoi

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
 * Texture view that is used for drawing the particles on canvas. The view is using
 * HA(hardware acceleration) for better performance.
 */
class TextureView : TextureView, SurfaceTextureListener {

    constructor(context: Context?) : super(context)
    constructor(context: Context?, attrs: AttributeSet?) : super(context, attrs)
    constructor(context: Context?, attrs: AttributeSet?, defStyleAttr: Int) : super(
        context,
        attrs,
        defStyleAttr
    )

    private var drawingThread: DrawingThread? = null
    lateinit var delaunay: Delaunay
    lateinit var voronoi: Voronoi
    lateinit var cells: Cells

    init {
        surfaceTextureListener = this
        isOpaque = false
    }

    /**
     * Create the drawing thread, attach the particles and start the thread.
     */
    override fun onSurfaceTextureAvailable(surface: SurfaceTexture, width: Int, height: Int) {

        initVoronoi(width, height)

        drawingThread = DrawingThread(
            TextureViewHolder(this),
            cells
        ).also {
            it.isRunning = true
            it.start()
        }

        this.setOnTouchListener(drawingThread)
    }

    /**
     * Update the view size in the particle properties.
     */
    override fun onSurfaceTextureSizeChanged(surface: SurfaceTexture, width: Int, height: Int) {

    }

    /**
     * Deestroy the drawing thread.
     */
    override fun onSurfaceTextureDestroyed(surface: SurfaceTexture): Boolean {
        var retry = true
        drawingThread?.isRunning = false
        while (retry) {
            try {
                drawingThread?.join()
                retry = false
            } catch (e: InterruptedException) {
                e.printStackTrace()
            }
        }
        return false
    }

    override fun onSurfaceTextureUpdated(surface: SurfaceTexture) {}

    private fun initVoronoi(width: Int, height: Int) {

        val numberOfRandomPoints = 50
        val relaxationLoops = 10

        // generate random points
        val points = ArrayList<PointF>()
        for (i in 0 until numberOfRandomPoints) {
            points.add(
                PointF(
                    (Math.random() * width - 1).toFloat(),
                    (Math.random() * height - 1).toFloat()
                )
            )
        }

        // generate delaunay and voronoi objects
        delaunay = Delaunay(points)
        voronoi = Voronoi(delaunay, RectF(0.0f, 0.0f, width.toFloat(), height.toFloat()))

        // apply relaxation to the points
        voronoi.relax(relaxationLoops)

        cells = Cells(
            doColorFilter(BitmapFactory.decodeResource(resources, R.drawable.cell3), 255, 255, 0, 120),
            14, 14,
            voronoi
        )
    }


    /**
     * Change the color of a bitmap by applying a color layer above this one, and
     * merging the pixels.
     */
    fun doColorFilter(src: Bitmap, red: Int, green: Int, blue: Int, alpha: Int): Bitmap {

        // image size
        val width = src.width
        val height = src.height

        // create output bitmap
        val bmOut = Bitmap.createBitmap(width, height, src.config)
        var pixel: Int

        // scan through all pixels
        for (x in 0 until width) {
            for (y in 0 until height) {

                // get pixel color
                pixel = src.getPixel(x, y)

                // apply filtering on each channel R, G, B
                val aB = Color.alpha(pixel)
                val rB = Color.red(pixel)
                val gB = Color.green(pixel)
                val bB = Color.blue(pixel)

                val aA = alpha
                val rA = red
                val gA = green
                val bA = blue

                val rOut = (rA * aA / 255) + (rB * aB * (255 - aA) / (255 * 255))
                val gOut = (gA * aA / 255) + (gB * aB * (255 - aA) / (255 * 255))
                val bOut = (bA * aA / 255) + (bB * aB * (255 - aA) / (255 * 255))
                var aOut = aA + (aB * (255 - aA) / 255)
                if (aB == 0) {
                    aOut = 0
                }

                // set new color pixel to output bitmap
                bmOut.setPixel(x, y, Color.argb(aOut, rOut, gOut, bOut))
            }
        }

        // return final image
        return bmOut
    }

}