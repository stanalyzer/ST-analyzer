/* ------------------------------------------------------
 * Function to update project grid by JAYSON
 * ------------------------------------------------------*/
function PrjUpdateGridByJson(obj) {
    var data = [];
    var indent = 0;
    var parents = [];
    var field_parent = [];
    var title = [];
    var path1 = [];
    var path2 = [];
    var path3 = [];
    var path4 = [];
    var path5 = [];
    var pbs = [];
    var date = [];
    var pkey = [];

    for (var i = 0; i < obj.parents.length; i++) {
      parents.push(obj.parents[i]);
    }


    for (var i = 0; i < obj.parent.length; i++) {
      field_parent.push(obj.parent[i]);
      title.push(obj.title[i]);
      path1.push(obj.path1[i]);
      path2.push(obj.path2[i]);
      path3.push(obj.path3[i]);
      path4.push(obj.path4[i]);
      path5.push(obj.path5[i]);
      pbs.push(obj.pbs[i]);
      date.push(obj.date[i]);
      pkey.push(obj.pkey[i]);
    }

    for (var i = 0; i < title.length; i++) {
      var d = (data[i] = {});
      var parent = null;
      
      if (field_parent[i] >= 0) {
        indent = 1;
        parent = field_parent[i];
      } else {
          parent = null;
          indent = 0;
      }
      d["id"] = "id_" + i;
      d["pkey"] = pkey[i];
      d["indent"] = indent;
      d["parent"] = parent;
      d["title"] = title[i];
      d["path1"] = path1[i];
      d["path2"] = path2[i];
      d["path3"] = path3[i];
      d["path4"] = path4[i];
      d["path5"] = path5[i];
      d["pbs"] = pbs[i];
      d["date"] = date[i];
    }
    
    dataView = new Slick.Data.DataView({ inlineFilters: true });
    
    grid = new Slick.Grid("#prjGrid", dataView, columns, options);
    
    grid.onClick.subscribe(function (e, args) {
      if ($j(e.target).hasClass("toggle")) {
        var item = dataView.getItem(args.row);
        if (item) {
          if (!item._collapsed) {
            item._collapsed = true;
          } else {
            item._collapsed = false;
          }
          dataView.updateItem(item.id, item);
        }
        e.stopImmediatePropagation();
      }
    });
    
    grid.setSelectionModel(new Slick.RowSelectionModel());
    
    var pager = new Slick.Controls.Pager(dataView, grid, $j("#pager"));
    var columnpicker = new Slick.Controls.ColumnPicker(columns, grid, options);
    
    
    // move the filter panel defined in a hidden div into grid top panel
    
    $j("#searchBox")
      .appendTo(grid.getTopPanel())
      .show();
    
    grid.onCellChange.subscribe(function (e, args) {
    dataView.updateItem(args.item.id, args.item);
    });
    
    
    grid.onAddNewRow.subscribe(function (e, args) {
    var item = {"num": data.length, "id": "new_" + (Math.round(Math.random() * 10000)), "title": "New task", "duration": "1 day", "percentComplete": 0, "start": "01/01/2009", "finish": "01/01/2009", "effortDriven": false};
    $j.extend(item, args.item);
    dataView.addItem(item);
    });
    
    grid.onKeyDown.subscribe(function (e) {
    // select all rows on ctrl-a
    if (e.which != 65 || !e.ctrlKey) {
      return false;
    }
    
    var rows = [];
    for (var i = 0; i < dataView.getLength(); i++) {
      rows.push(i);
    }
    
    grid.setSelectedRows(rows);
    e.preventDefault();
    });
    
    grid.onSort.subscribe(function (e, args) {
    sortdir = args.sortAsc ? 1 : -1;
    sortcol = args.sortCol.field;
    
    if ($j.browser.msie && $j.browser.version <= 8) {
      // using temporary Object.prototype.toString override
      // more limited and does lexicographic sort only by default, but can be much faster
      // use numeric sort of % and lexicographic for everything else
      dataView.fastSort(sortcol, args.sortAsc);
    } else {
      // using native sort with comparer
      // preferred method but can be very slow in IE with huge datasets
      dataView.sort(comparer, args.sortAsc);
    }
    });
    
    // wire up model events to drive the grid
    dataView.onRowCountChanged.subscribe(function (e, args) {
    grid.updateRowCount();
    grid.render();
    });
    
    dataView.onRowsChanged.subscribe(function (e, args) {
    grid.invalidateRows(args.rows);
    grid.render();
    });
    
    dataView.onPagingInfoChanged.subscribe(function (e, pagingInfo) {
    var isLastPage = pagingInfo.pageNum == pagingInfo.totalPages - 1;
    var enableAddRow = isLastPage || pagingInfo.pageSize == 0;
    var options = grid.getOptions();
    
    if (options.enableAddRow != enableAddRow) {
      grid.setOptions({enableAddRow: enableAddRow});
    }
    });
    
    
    var h_runfilters = null;
    
    // wire up the search textbox to apply the filter to the model
    $j("#txtSearch").keyup(function (e) {
    Slick.GlobalEditorLock.cancelCurrentEdit();
    
    // clear on Esc
    if (e.which == 27) {
      this.value = "";
    }
    
    searchString = this.value;
    updateFilter();
    });
    
    function updateFilter() {
    dataView.setFilterArgs({
      searchString: searchString
    });
    dataView.refresh();
    }
    
    // initialize the model after all the events have been hooked up
    dataView.beginUpdate();
    dataView.setItems(data);
    dataView.setFilterArgs({
        searchString: searchString
    });
    dataView.setFilter(myFilter);
    dataView.endUpdate();
    
    
    // if you don't want the items that are not visible (due to being filtered out
    // or being on a different page) to stay selected, pass 'false' to the second arg
    dataView.syncGridSelection(grid, true);
}



/* ------------------------------------------------------
 * Function to update user grid by JAYSON
 * ------------------------------------------------------*/
function UsrUpdateGridByJson(obj) {
    var data = [];
    var indent = 0;
    var parents = [];
    var field_parent = [];
    var uid = [];
    var pwd = [];
    var email = [];
    var level = [];
    var pkey = [];
    
    for (var i = 0; i < obj.parents.length; i++) {
      parents.push(obj.parents[i]);
    }

    for (var i = 0; i < obj.parent.length; i++) {
      field_parent.push(obj.parent[i]);
      uid.push(obj.uid[i]);
      pwd.push(obj.pwd[i]);
      email.push(obj.email[i]);
      level.push(obj.level[i]);
      pkey.push(obj.pkey[i]);
    }

    for (var i = 0; i < uid.length; i++) {
      var d = (data[i] = {});
      var parent = null;
      
      if (field_parent[i] >= 0) {
        indent = 1;
        parent = field_parent[i];
      } else {
          parent = null;
          indent = 0;
      }
      d["id"] = "id_" + i;
      d["pkey"] = pkey[i];
      d["indent"] = indent;
      d["parent"] = parent;
      d["uid"] = uid[i];
      d["pwd"] = pwd[i];
      d["email"] = email[i];
      d["level"] = level[i];
    }
    
    dataView = new Slick.Data.DataView({ inlineFilters: true });
    
    grid = new Slick.Grid("#usrGrid", dataView, columns, options);
    
    grid.onClick.subscribe(function (e, args) {
      if ($j(e.target).hasClass("toggle")) {
        var item = dataView.getItem(args.row);
        if (item) {
          if (!item._collapsed) {
            item._collapsed = true;
          } else {
            item._collapsed = false;
          }
          dataView.updateItem(item.id, item);
        }
        e.stopImmediatePropagation();
      }
    });
    
    grid.setSelectionModel(new Slick.RowSelectionModel());
    
    var pager = new Slick.Controls.Pager(dataView, grid, $j("#pager"));
    var columnpicker = new Slick.Controls.ColumnPicker(columns, grid, options);
    
    
    // move the filter panel defined in a hidden div into grid top panel
    
    $j("#searchBox")
      .appendTo(grid.getTopPanel())
      .show();
    
    grid.onCellChange.subscribe(function (e, args) {
    dataView.updateItem(args.item.id, args.item);
    });
    
    
    grid.onAddNewRow.subscribe(function (e, args) {
    var item = {"num": data.length, "id": "new_" + (Math.round(Math.random() * 10000)), "title": "New task", "duration": "1 day", "percentComplete": 0, "start": "01/01/2009", "finish": "01/01/2009", "effortDriven": false};
    $j.extend(item, args.item);
    dataView.addItem(item);
    });
    
    grid.onKeyDown.subscribe(function (e) {
    // select all rows on ctrl-a
    if (e.which != 65 || !e.ctrlKey) {
      return false;
    }
    
    var rows = [];
    for (var i = 0; i < dataView.getLength(); i++) {
      rows.push(i);
    }
    
    grid.setSelectedRows(rows);
    e.preventDefault();
    });
    
    grid.onSort.subscribe(function (e, args) {
    sortdir = args.sortAsc ? 1 : -1;
    sortcol = args.sortCol.field;
    
    if ($j.browser.msie && $j.browser.version <= 8) {
      // using temporary Object.prototype.toString override
      // more limited and does lexicographic sort only by default, but can be much faster
      // use numeric sort of % and lexicographic for everything else
      dataView.fastSort(sortcol, args.sortAsc);
    } else {
      // using native sort with comparer
      // preferred method but can be very slow in IE with huge datasets
      dataView.sort(comparer, args.sortAsc);
    }
    });
    
    // wire up model events to drive the grid
    dataView.onRowCountChanged.subscribe(function (e, args) {
    grid.updateRowCount();
    grid.render();
    });
    
    dataView.onRowsChanged.subscribe(function (e, args) {
    grid.invalidateRows(args.rows);
    grid.render();
    });
    
    dataView.onPagingInfoChanged.subscribe(function (e, pagingInfo) {
    var isLastPage = pagingInfo.pageNum == pagingInfo.totalPages - 1;
    var enableAddRow = isLastPage || pagingInfo.pageSize == 0;
    var options = grid.getOptions();
    
    if (options.enableAddRow != enableAddRow) {
      grid.setOptions({enableAddRow: enableAddRow});
    }
    });
    
    
    var h_runfilters = null;
    
    // wire up the search textbox to apply the filter to the model
    $j("#txtSearch").keyup(function (e) {
    Slick.GlobalEditorLock.cancelCurrentEdit();
    
    // clear on Esc
    if (e.which == 27) {
      this.value = "";
    }
    
    searchString = this.value;
    updateFilter();
    });
    
    function updateFilter() {
    dataView.setFilterArgs({
      searchString: searchString
    });
    dataView.refresh();
    }
    
    // initialize the model after all the events have been hooked up
    dataView.beginUpdate();
    dataView.setItems(data);
    dataView.setFilterArgs({
        searchString: searchString
    });
    dataView.setFilter(myFilter);
    dataView.endUpdate();
    
    
    // if you don't want the items that are not visible (due to being filtered out
    // or being on a different page) to stay selected, pass 'false' to the second arg
    dataView.syncGridSelection(grid, true);
}
