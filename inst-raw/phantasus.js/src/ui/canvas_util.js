phantasus.CanvasUtil = function () {
};
phantasus.CanvasUtil.dragging = false;

phantasus.CanvasUtil.FONT_NAME = '"Helvetica Neue",Helvetica,Arial,sans-serif';
phantasus.CanvasUtil.FONT_COLOR = 'rgb(34, 34, 34)';
phantasus.CanvasUtil.getPreferredSize = function (c) {
  var size = c.getPreferredSize();
  var prefWidth = c.getPrefWidth();
  var prefHeight = c.getPrefHeight();
  // check for override override
  if (prefWidth !== undefined) {
    size.width = prefWidth;
  }
  if (prefHeight !== undefined) {
    size.height = prefHeight;
  }
  return size;
};
phantasus.CanvasUtil.BACKING_SCALE = 1;
if (typeof window !== 'undefined' && 'devicePixelRatio' in window) {
  if (window.devicePixelRatio > 1) {
    phantasus.CanvasUtil.BACKING_SCALE = window.devicePixelRatio;
  }
}

phantasus.CanvasUtil.setBounds = function (canvas, bounds) {
  var backingScale = phantasus.CanvasUtil.BACKING_SCALE;

  if (bounds.height != null) {
    canvas.height = bounds.height * backingScale;
    canvas.style.height = bounds.height + 'px';
  }
  if (bounds.width != null) {
    canvas.width = bounds.width * backingScale;
    canvas.style.width = bounds.width + 'px';
  }
  if (bounds.left != null) {
    canvas.style.left = bounds.left + 'px';
  }
  if (bounds.top != null) {
    canvas.style.top = bounds.top + 'px';
  }
};

phantasus.CanvasUtil.drawShape = function (context, shape, x, y, size2) {
  if (size2 < 0) {
    return;
  }
  context.beginPath();
  if (shape === 'minus') {
    context.arc(x, y, size2, 0, 2 * Math.PI, false);
    context.moveTo(x - size2, y);
    context.lineTo(x + size2, y);
  } else if (shape === 'circle') {
    context.arc(x, y, size2, 0, 2 * Math.PI, false);
  } else if (shape === 'square') {
    context.rect(x - size2, y - size2, size2 * 2, size2 * 2);
  } else if (shape === 'plus') {
    // vertical line
    context.moveTo(x, y - size2);
    context.lineTo(x, y + size2);
    // horizontal line
    context.moveTo(x - size2, y);
    context.lineTo(x + size2, y);
  } else if (shape === 'x') {
    context.moveTo(x - size2, y - size2);
    context.lineTo(x + size2, y + size2);
    context.moveTo(x + size2, y - size2);
    context.lineTo(x - size2, y + size2);
  } else if (shape === 'asterisk') {
    // x with vertical line
    context.moveTo(x - size2, y - size2);
    context.lineTo(x + size2, y + size2);
    context.moveTo(x + size2, y - size2);
    context.lineTo(x - size2, y + size2);

    context.moveTo(x, y - size2);
    context.lineTo(x, y + size2);
  } else if (shape === 'diamond') {
    // start at middle top
    context.moveTo(x, y - size2);
    // right
    context.lineTo(x + size2, y);
    // bottom
    context.lineTo(x, y + size2);
    // left
    context.lineTo(x - size2, y);
    // top
    context.lineTo(x, y - size2);
  } else if (shape === 'triangle-up') {
    // top
    context.moveTo(x, y - size2);
    // right
    context.lineTo(x + size2, y + size2);
    // left
    context.lineTo(x - size2, y + size2);
    context.lineTo(x, y - size2);
  } else if (shape === 'triangle-down') {
    // bottom
    context.moveTo(x, y + size2);
    // left
    context.lineTo(x - size2, y - size2);
    // right
    context.lineTo(x + size2, y - size2);
    context.lineTo(x, y + size2);
  } else if (shape === 'triangle-left') {
    // left
    context.moveTo(x - size2, y);
    // top
    context.lineTo(x + size2, y - size2);
    // bottom
    context.lineTo(x + size2, y + size2);
    context.lineTo(x - size2, y);
  } else if (shape === 'triangle-right') {
    // right
    context.moveTo(x + size2, y);
    // lower left
    context.lineTo(x - size2, y + size2);

    // upper left
    context.lineTo(x - size2, y - size2);
    context.lineTo(x + size2, y);
  }
  context.stroke();

};
phantasus.CanvasUtil.drawLine = function (context, x1, y1, x2, y2) {
  context.beginPath();
  context.moveTo(x1, y1);
  context.lineTo(x2, y2);
  context.stroke();
};
phantasus.CanvasUtil.resetTransform = function (context) {
  context.setTransform(1, 0, 0, 1, 0, 0);
  if (phantasus.CanvasUtil.BACKING_SCALE !== 1) {
    context.scale(phantasus.CanvasUtil.BACKING_SCALE,
      phantasus.CanvasUtil.BACKING_SCALE);
  }
};
phantasus.CanvasUtil.bezierCurveTo = function (context, start, end) {
  var m1 = (start[1] + end[1]) / 2;
  context.beginPath();
  context.moveTo(start[0], start[1]);
  // context.lineTo(leftp[0], leftp[1]);
  context.bezierCurveTo(start[0], m1, end[0], m1, end[0], end[1]);
  context.stroke();
};
phantasus.CanvasUtil.createCanvas = function () {
  var $c = $('<canvas></canvas>');
  $c.attr('tabindex', '0');
  $c.css({
    cursor: 'default',
    outline: 0,
    overflow: 'hidden',
    position: 'absolute',
    'z-index': 1
  });
  return $c[0];
};
phantasus.CanvasUtil.getHeaderStringWidth = function (context, s) {
  context.font = '14px ' + phantasus.CanvasUtil.FONT_NAME;
  return context.measureText(s).width + 18;
};
phantasus.CanvasUtil.getVectorStringWidth = function (context, vector, positions,
                                                     end) {
  if (positions.getSize() < 6) {
    return 0;
  }
  var fontSize = Math.min(24, positions.getSize() - 2);
  if (fontSize <= 0) {
    return 0;
  }

  context.font = fontSize + 'px ' + phantasus.CanvasUtil.FONT_NAME;

  var toString = phantasus.VectorTrack.vectorToString(vector);
  var maxWidth = 0;
  // var maxWidth2 = 0;
  var n = end <= 0 ? vector.size() : Math.min(end, vector.size());
  for (var i = 0; i < n; i++) {
    var value = vector.getValue(i);
    if (value != null && value != '') {
      value = toString(value);
    } else {
      continue;
    }
    var width = context.measureText(value).width;
    if (width > maxWidth) {
      maxWidth = width;
    }
    // if (width > maxWidth2 && width < maxWidth) {
    // maxWidth2 = width;
    // }
  }
  return maxWidth > 0 ? (maxWidth + 2) : maxWidth;
};
phantasus.CanvasUtil.clipString = function (context, string, availTextWidth) {
  var textWidth = context.measureText(string).width;
  if (textWidth <= availTextWidth) {
    return string;
  }
  var clipString = '...';
  availTextWidth -= context.measureText(clipString).width;
  if (availTextWidth <= 0) {
    // can not fit any characters
    return clipString;
  }
  var width = 0;
  for (var nChars = 0, stringLength = string.length; nChars < stringLength; nChars++) {
    width += context.measureText(string[nChars]).width;
    if (width > availTextWidth) {
      string = string.substring(0, nChars);
      break;
    }
  }
  return string + clipString;
};
phantasus.CanvasUtil.toSVG = function (drawable, file) {
  var totalSize = {
    width: drawable.getWidth(),
    height: drawable.getHeight()
  };
  var context = new C2S(totalSize.width, totalSize.height);
  context.save();
  drawable.draw({
    x: 0,
    y: 0,
    width: totalSize.width,
    height: totalSize.height
  }, context);
  context.restore();
  var svg = context.getSerializedSvg();
  var blob = new Blob([svg], {
    type: 'text/plain;charset=utf-8'
  });
  saveAs(blob, file);
};
phantasus.CanvasUtil.getMousePos = function (element, event, useDelta) {
  return phantasus.CanvasUtil.getMousePosWithScroll(element, event, 0, 0,
    useDelta);
};

phantasus.CanvasUtil.getClientXY = function (event, useDelta) {
  var clientX;
  var clientY;
  if (event.pointers) {
    if (event.pointers.length > 0) {
      clientX = event.pointers[0].clientX - (useDelta ? event.deltaX : 0);
      clientY = event.pointers[0].clientY - (useDelta ? event.deltaY : 0);
    } else {
      clientX = event.srcEvent.clientX - (useDelta ? event.deltaX : 0);
      clientY = event.srcEvent.clientY - (useDelta ? event.deltaY : 0);
    }
  } else {
    clientX = event.clientX;
    clientY = event.clientY;
  }
  return {
    x: clientX,
    y: clientY
  };
};
phantasus.CanvasUtil.getMousePosWithScroll = function (element, event, scrollX,
                                                      scrollY, useDelta) {
  return phantasus.CanvasUtil._getMousePosWithScroll(element, scrollX,
    scrollY, phantasus.CanvasUtil.getClientXY(event, useDelta));
};

phantasus.CanvasUtil._getMousePosWithScroll = function (element, scrollX,
                                                       scrollY, clientXY) {
  var rect = element.getBoundingClientRect();
  return {
    x: clientXY.x - rect.left + scrollX,
    y: clientXY.y - rect.top + scrollY
  };
};
