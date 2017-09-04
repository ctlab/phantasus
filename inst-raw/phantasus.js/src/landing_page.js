phantasus.LandingPage = function (pageOptions) {
  pageOptions = $.extend({}, {
    el: $('#vis')
  }, pageOptions);
  this.pageOptions = pageOptions;
  var _this = this;

  var $el = $('<div class="container" style="display: none;"></div>');
  this.$el = $el;
  var html = [];
  html.push('<div data-name="help" class="pull-right"></div>');

  html
    .push('<div style="margin-bottom:10px;"><svg width="32px" height="32px"><g><rect x="0" y="0" width="32" height="14" style="fill:#ca0020;stroke:none"/><rect x="0" y="18" width="32" height="14" style="fill:#0571b0;stroke:none"/></g></svg> <div data-name="brand" style="display:inline-block; vertical-align: top;font-size:24px;font-family:sans-serif;">');
  html.push('<span>P</span>');
  html.push('<span>h</span>');
  html.push('<span>a</span>');
  html.push('<span>n</span>');
  html.push('<span>t</span>');
  html.push('<span>a</span>');
  html.push('<span>s</span>');
  html.push('<span>u</span>');
  html.push('<span>s</span>');
  html.push('</span>');
  html.push('</div>');

  html.push('<h4>Open your own file</h4>');
  html.push('<div data-name="formRow" class="center-block"></div>');
  html.push('<div style="display: none;" data-name="preloadedDataset"><h4>Or select a preloaded' +
    ' dataset</h4></div>');
  html.push('</div>');
  var $html = $(html.join(''));
  var colorScale = d3.scale.linear().domain([0, 4, 7]).range(['#ca0020', '#999999', '#0571b0']).clamp(true);
  var brands = $html.find('[data-name="brand"] > span');
  $html.appendTo($el);
  new phantasus.HelpMenu().$el.appendTo($el.find('[data-name=help]'));
  var formBuilder = new phantasus.FormBuilder();
  formBuilder.append({
    name: 'file',
    showLabel: false,
    value: '',
    type: 'file',
    required: true,
    help: phantasus.DatasetUtil.DATASET_FILE_FORMATS
  });
  formBuilder.$form.appendTo($el.find('[data-name=formRow]'));
  this.formBuilder = formBuilder;
  this.$sampleDatasetsEl = $el.find('[data-name=preloadedDataset]');
  var index = 0;
  var step = function () {
    brands[index].style.color = colorScale(index);
    index++;
    if (index < brands.length) {
      setTimeout(step, 200);
    }
  };
  setTimeout(step, 300);
  this.tabManager = new phantasus.TabManager({landingPage: this});
  this.tabManager.on('change rename add remove', function (e) {
    var title = _this.tabManager.getTabText(_this.tabManager.getActiveTabId());
    if (title == null || title === '') {
      title = 'phantasus';
    }
    document.title = title;
  });

  this.tabManager.$nav.appendTo($(this.pageOptions.el));
  this.tabManager.$tabContent.appendTo($(this.pageOptions.el));
  // for (var i = 0; i < brands.length; i++) {
  // 	brands[i].style.color = colorScale(i);
  // }
};

phantasus.LandingPage.prototype = {
  open: function (openOptions) {
    this.dispose();
    var optionsArray = _.isArray(openOptions) ? openOptions : [openOptions];
    var _this = this;
    // console.log(optionsArray);
    for (var i = 0; i < optionsArray.length; i++) {
      var options = optionsArray[i];
      options.tabManager = _this.tabManager;
      options.focus = i === 0;
      options.landingPage = _this;

      if (options.dataset.options.isGEO) {
        var req = ocpu.call('checkGPLs', { name : options.dataset.file }, function (session) {
          session.getObject(function (filenames) {
            filenames = JSON.parse(filenames);
            // console.log(filenames);
            if (filenames.length === 0) {
              alert("Dataset" + " " + options.dataset.file + " does not exist");
              _this.show();
            }
            if (filenames.length === 1) {
              new phantasus.HeatMap(options);
            }
            else {
              for (var j = 0; j < filenames.length; j++) {
                var specificOptions = options;
                specificOptions.dataset.file = filenames[j];

                new phantasus.HeatMap(specificOptions);
              }
            }
          })
        });
        req.fail(function () {
          new Error("Checking GPLs call to OpenCPU failed" + req.responseText);
        });
      }
      else {
        new phantasus.HeatMap(options);
      }
    }

  },
  dispose: function () {
    this.formBuilder.setValue('file', '');
    this.$el.hide();
    $(window)
      .off(
        'paste.phantasus drop.phantasus dragover.phantasus dragenter.phantasus');
    this.formBuilder.off('change');
  },
  show: function () {
    var _this = this;
    if (navigator.onLine && !this.sampleDatasets) {
      this.sampleDatasets = new phantasus.SampleDatasets({
        $el: this.$sampleDatasetsEl,
        show: true,
        callback: function (heatMapOptions) {
          _this.open(heatMapOptions);
        }
      });
    }

    this.$el.show();

    this.formBuilder.on('change', function (e) {
      var value = e.value;
      if (value !== '' && value != null) {
        _this.openFile(value);
      }
    });

    $(window).on('beforeunload.phantasus', function () {
      if (_this.tabManager.getTabCount() > 0) {
        return 'Are you sure you want to close phantasus?';
      }
    });
    $(window).on('paste.phantasus', function (e) {
      var tagName = e.target.tagName;
      if (tagName == 'INPUT' || tagName == 'SELECT' || tagName == 'TEXTAREA') {
        return;
      }

      var text = e.originalEvent.clipboardData.getData('text/plain');
      if (text != null && text.length > 0) {
        e.preventDefault();
        e.stopPropagation();
        var url;
        if (text.indexOf('http') === 0) {
          url = text;
        } else {
          var blob = new Blob([text]);
          url = window.URL.createObjectURL(blob);
        }

        _this.openFile(url);
      }

    }).on('dragover.phantasus dragenter.phantasus', function (e) {
      e.preventDefault();
      e.stopPropagation();
    }).on(
      'drop.phantasus',
      function (e) {
        if (e.originalEvent.dataTransfer
          && e.originalEvent.dataTransfer.files.length) {
          e.preventDefault();
          e.stopPropagation();
          var files = e.originalEvent.dataTransfer.files;
          _this.openFile(files[0]);
        } else if (e.originalEvent.dataTransfer) {
          var url = e.originalEvent.dataTransfer.getData('URL');
          e.preventDefault();
          e.stopPropagation();
          _this.openFile(url);
        }
      });
  },
  openFile: function (value) {
    var _this = this;
    var fileName = phantasus.Util.getFileName(value);
    if (fileName.toLowerCase().indexOf('.json') === fileName.length - 5) {
      phantasus.Util.getText(value).done(function (text) {
        _this.open(JSON.parse(text));
      }).fail(function (err) {
        phantasus.FormBuilder.showMessageModal({
          title: 'Error',
          message: 'Unable to load session'
        });
      });
    } else {
      var options = {
        dataset: {
          file: value,
          options: {
            interactive: true,
            isGEO: (fileName.toUpperCase().indexOf('GSE') === 0 || fileName.toUpperCase().indexOf('GDS') === 0) && fileName.indexOf('.') === -1
          }
        }
      };

      phantasus.OpenDatasetTool.fileExtensionPrompt(fileName, function (readOptions) {
        if (readOptions) {
          for (var key in readOptions) {
            options.dataset.options[key] = readOptions[key];
          }
        }
        _this.open(options);
      });
    }
  }
};
