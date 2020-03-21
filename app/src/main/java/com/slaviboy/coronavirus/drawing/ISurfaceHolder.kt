package com.slaviboy.coronavirus.drawing
import android.graphics.Canvas

interface ISurfaceHolder {
    fun unlockCanvasAndPost(canvas: Canvas?)
    fun lockCanvas(): Canvas?
}