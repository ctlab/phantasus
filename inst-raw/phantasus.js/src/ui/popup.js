phantasus.Popup = {};
phantasus.Popup.initted = false;
phantasus.Popup.init = function () {
  if (phantasus.Popup.initted) {
    return;
  }

  phantasus.Popup.initted = true;
  phantasus.Popup.$popupDiv = $(document.createElement('div'));
  phantasus.Popup.$popupDiv.css('position', 'absolute').css('zIndex', 1050).css('overflow', 'auto').addClass('dropdown clearfix');
  phantasus.Popup.$contextMenu = $(document.createElement('ul'));
  phantasus.Popup.$contextMenu.addClass('dropdown-menu').css('display',
    'block').css('position', 'static').css('margin-bottom', '5px');
  phantasus.Popup.$contextMenu.appendTo(phantasus.Popup.$popupDiv);
  phantasus.Popup.$contextMenu.on('click', 'a', function (e) {
    e.preventDefault();
    var $this = $(this);
    // if (!$this.hasClass('copy')) {
    phantasus.Popup.popupCallback(e, $this.data('name'));
    phantasus.Popup.hide();
    // }

  });
};

phantasus.Popup.popupInDom = false;
phantasus.Popup.hidePopupMenu = function (e) {
  if (phantasus.Popup.component == e.target) {
    e.preventDefault();
    e.stopPropagation();
  }
  phantasus.Popup.hide();
};
phantasus.Popup.hide = function () {
  phantasus.Popup.$popupDiv.hide();
  $(document.body).off('mousedown', phantasus.Popup.hidePopupMenu);
  phantasus.Popup.popupCallback = null;
  phantasus.Popup.component = null;
};

phantasus.Popup.showPopup = function (menuItems, position, component, callback) {
  phantasus.Popup.init();
  if (phantasus.Popup.component == component) {
    phantasus.Popup.hide();
    return;
  }
  phantasus.Popup.popupCallback = callback;
  phantasus.Popup.component = component;
  var html = [];
  for (var i = 0, length = menuItems.length; i < length; i++) {
    var item = menuItems[i];
    if (item.header) {
      html.push('<li role="presentation" class="dropdown-header">'
        + item.name + '</li>');
    } else if (item.separator) {
      html.push('<li class="divider"></li>');
    } else {
      html.push('<li role="presentation"');
      if (item.disabled) {
        html.push('class="disabled"');
      }
      html.push('><a data-name="' + item.name
        + '" data-type="popup-item" tabindex="-1" href="#"');
      if (item.class) {
        html.push(' class="' + item.class + '"');
      }
      html.push('>');
      if (item.checked) {
        html
          .push('<span class="dropdown-checkbox fa fa-check"></span>');
      }

      html.push(item.name);
      if (item.icon) {
        html.push('<span class="pull-right ' + item.icon + '"></span>');
      }
      html.push('</a>');

      html.push('</li>');
    }
  }
  phantasus.Popup.$contextMenu.html(html.join(''));
  if (!phantasus.Popup.popupInDom) {
    phantasus.Popup.popupInDom = true;
    phantasus.Popup.$popupDiv.appendTo($(document.body));
  }
  var $body = $(document.body);
  var $window = $(window);
  var windowWidth = $window.width();
  var windowHeight = $window.height();
  var popupWidth = phantasus.Popup.$popupDiv.width();
  var popupHeight = phantasus.Popup.$popupDiv.height();
  var left = position.x;
  var top = position.y;
  // default is bottom-right
  if ((left + popupWidth) >= windowWidth) { // offscreen right
    left -= popupWidth;
    left = Math.max(4, left);
  }
  if ((top + popupHeight) >= (windowHeight)) { // offscreen bottom
    top -= popupHeight;
    top = Math.max(4, top);
  }

  phantasus.Popup.$popupDiv.css({
    height: popupHeight + 'px',
    display: 'block',
    left: left,
    top: top
  });

  phantasus.Popup.$popupDiv.show();

  $body.off('mousedown', phantasus.Popup.hidePopupMenu);
  window.setTimeout(function () {
    $body.on('mousedown', function (e) {
      var $target = $(e.target);
      if ($target[0] !== phantasus.Popup.$popupDiv[0] && $target.data('type') !== 'popup-item') {
        phantasus.Popup.hidePopupMenu(e);
      }
    });
  }, 1);
};
