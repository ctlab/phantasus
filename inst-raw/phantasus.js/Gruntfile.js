module.exports = function (grunt) {
  grunt
  .initConfig({
    pkg: grunt.file.readJSON('package.json'),
    meta: {
      banner: '/*! <%= pkg.title || pkg.name %> - v<%= pkg.version %> - '
      + '<%= grunt.template.today("yyyy-mm-dd") %>\n'
      + '<%= pkg.homepage ? "* " + pkg.homepage + "\n" : "" %>'
      + '* Copyright (c) <%= grunt.template.today("yyyy") %> <%= pkg.author.name %>; %> */'
    },
    uglify: {
      phantasus: {
        options: {
          mangle: false,
          compress: false,
          preserveComments: false
        },
        files: {
          'js/phantasus-latest.min.js': ['js/phantasus.js']
        }
      },
      extJs: {
        options: {
          mangle: false,
          compress: false,
          preserveComments: false
        },
        files: {
          'js/phantasus-external.min.js': ['js/phantasus-external.js']
        }
      }
    },
    cssmin: {
      css: {
        src: 'css/phantasus.all.css',
        dest: 'css/phantasus-latest.min.css'
      }
    },
    concat: {
      css: {
        src: ['css/bootstrap.min.css',
          'css/bootstrap-select.min.css',
          'css/jquery-ui.min.css',
          'css/font-awesome.min.css',
          'css/slick.grid.css', 'css/phantasus.grid.css',
          'css/animate.css', 'css/phantasus.css'],
        dest: 'css/phantasus.all.css'
      },
      extJsAll: {
        src: ['js/phantasus-external.min.js',
          'js/plotly-latest.min.js', 'js/papaparse.min.js'],
        dest: 'js/phantasus-external-latest.min.js'
      },
      extJs: {
        dest: 'js/phantasus-external.js',
        src: ['js/d3.min.js', 'js/jquery-2.2.4.min.js',
          'js/bootstrap.min.js', 'js/underscore-min.js',
          'js/newick.js', 'js/hammer.min.js',
          'js/jquery.mousewheel.min.js',
          'js/bootstrap-select.min.js',
          'js/xlsx.full.min.js', 'js/canvas2svg.js',
          'js/canvg.js', 'js/rgbcolor.js',
          'js/jquery-ui.min.js', 'js/parser.js',
          'js/colorbrewer.js',
          'js/jquery.event.drag-2.2.js',
          'js/clipboard.min.js', 'js/slick.min.js', 'js/canvas-toBlob.js',
          'js/js.cookie.js','js/long.js', 'js/bytebuffer.js', 'js/protobuf.js',
          'js/opencpu-0.5.js']
      },
      phantasus: {
        options: {
          banner: '(function(global){\n\'use strict\';\n',
          footer: '\n})(typeof window !== \'undefined\' ? window : this);\n'
        },
        dest: 'js/phantasus.js',
        src: ['src/util/util.js', 'src/util/*.js',
          'src/io/*.js', 'src/matrix/vector_adapter.js',
          'src/matrix/*.js', 'src/*.js',
          'src/tools/*.js', 'src/ui/*.js', 'src/**/*.js', 'js/tsne.js']
      }
    },
    watch: {
      files: ['src/*.js', 'src/**/*.js'],
      tasks: ['concat:phantasus']
    }
  });

  // rebuild js and css:
  // grunt concat:phantasus concat:extJs uglify concat:extJsAll
  grunt.registerTask('default', 'watch');
  grunt.loadNpmTasks('grunt-contrib-watch');
  grunt.loadNpmTasks('grunt-contrib-concat');
  grunt.loadNpmTasks('grunt-contrib-uglify');
  grunt.loadNpmTasks('grunt-contrib-copy');
  grunt.loadNpmTasks('grunt-contrib-cssmin');
  grunt.loadNpmTasks('grunt-contrib-jshint');
};
