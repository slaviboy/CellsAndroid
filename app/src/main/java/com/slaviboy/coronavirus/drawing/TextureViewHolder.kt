package com.slaviboy.coronavirus.drawing

import android.graphics.Canvas

class TextureViewHolder(private val textureView: TextureView) :
    ISurfaceHolder {

    override fun unlockCanvasAndPost(canvas: Canvas?) {
        textureView.unlockCanvasAndPost(canvas)
    }

    override fun lockCanvas(): Canvas? {
        return textureView.lockCanvas()
    }

}