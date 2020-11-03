/*
* @Author: UnsignedByte
* @Date:   14:51:38, 09-Jun-2020
* @Last Modified by:   UnsignedByte
* @Last Modified time: 18:12:40, 02-Nov-2020
*/

data = $$datainput$$;

// data = JSON.parse(data);

//---------------------------------------------------- USED FOR GENERATING SHEETS FROM DATA ----------------------------------------------------

//generate a set of tables to simulate excel with sheets
function generateSheet(datas, names, types){
  footerElems.isColoured = false;
  let sheetbar = $("#sheetbar").empty();
  let tables = $("#tables").empty();
  for(let i = 0; i < datas.length; i++){
    $("<input/>", {type:"button", value:names[i]}).click(
      function(){
        footerElems.isColoured = false;
        $("#sheetbar input.selected").removeClass('selected');
        $(this).addClass('selected');
        $("#tables div.container.selected").removeClass('selected');
        $("#tables div.container").eq($(this).index()).addClass('selected');
      }).appendTo(sheetbar);

    let cont = $("<div/>", {class:"container clusterize"}).appendTo(tables);
    switch (types[i]){
      case 'basic':
        basicTableCreate(datas[i]).appendTo(cont);
        break;
      case 'class':
        classTableCreate(datas[i]).appendTo(cont);
        break;
      case 'html': //Raw html data
      	$(datas[i]).appendTo(cont);
      	break;
      default:
        break;
    }
  }
  sheetbar.children().first().addClass('selected');
  tables.children().first().addClass('selected');
}

//basic table
function basicTableCreate(dat){
  let table = $('<table/>', {class:'basicTable'});
  headers = Object.keys(dat[0]);
  $('<tr/>').addRow('th', [...Array(headers.length).keys()], ["Row Numbers", ...[...Array(headers.length).keys()].slice(1)]).wrap('<thead/>').parent().appendTo(table);
  $('<tr/>').addRow('th', [...Array(headers.length).keys()], headers).wrap('<thead/>').parent().appendTo(table);

  let body = $('<tbody/>').appendTo(table);
  for(let i = 0; i < dat.length; i++){
    if (i%1000 === 0) console.log(`${i} rows rendered`);
    $('<tr/>').addRow('td', headers, dat[i]).appendTo(body);
  }

  return table;
}


//table with classes (expandable)
function classTableCreate(dat){
  let table = $('<table/>', {class:'classTable'});

  table.on('click', 'tr.trclass .fa-chevron-down', function(){
    $(this).closest('tbody').toggleClass('tropen');
  });

  headers = Object.keys(dat[0]).slice(0,-1);
  $('<tr/>').addRow('th', [...Array(headers.length).keys()], ["Row Numbers", ...[...Array(headers.length).keys()].slice(1)])
        .prepend($('<th/>')).wrap('<thead/>').parent().appendTo(table);
  $('<tr/>').addRow('th', [...Array(headers.length).keys()], headers).prepend($('<th/>')).wrap('<thead/>').parent().appendTo(table);

  let headerRow = null;
  for(let i = 0; i < dat.length; i++){
    if (i%1000 === 0) console.log(`${i} rows rendered`);
    // console.log(dat[i].Row_Type);
    if (dat[i].Row_Type){ //is a class
      headerRow = $('<tr/>', {class:'trclass'}).addRow('td', headers, dat[i]).prepend($('<td><i class="fa fa-chevron-down"></i></td>')).wrap('<tbody/>').parent().appendTo(table);
      if (table.children().length%2 == 1) headerRow.addClass('odd');
    }else{
      $('<tr/>', {class:'trchild'}).addRow('td', headers, dat[i]).prepend($('<td/>')).appendTo(headerRow);
    }
  }
  return table;
}

//populate table row from data
$.fn.extend({addRow: function(type, headers, dat) {
  for(let i = 0; i < headers.length; i++){
    $(`<${type}/>`, {text:dat[headers[i]]}).appendTo(this);
  }
  return this;
}});

// function addRow(row, type, headers, dat){
//   for(let i = 0; i < headers.length; i++){
//     $(`<${type}/>`, {text:dat[headers[i]]}).appendTo(row);
//   }
//   return row
// }

function pickGradient(color1, color2, weight) {
    var w1 = weight;
    var w2 = 1 - w1;

    return (Math.round((color1>>16) * w1 + (color2>>16) * w2)<<16)
        + (Math.round((color1>>8)%256 * w1 + (color2>>8)%256 * w2)<<8)
        + Math.round(color1%256 * w1 + color2%256 * w2);
    // return rghex(rgb);
}

function rghex(col){
  return (col[0]<<16)+(col[1]<<8)+col[2];
}

// $('#tables').append(classTableCreate(data.HeatMap[0].Data));
// $('#tables').append(basicTableCreate(data.ModMapper[0].Sheets[0]));

  // <input type="radio" name="tabletype" value="HeatMap"> HeatMap
  // <input type="radio" name="tabletype" value="Modmapper"> Modmapper

//buttons for time intensive operations on click (will show processing bar)
function longButton(name, f){
  return $('<input/>', {type:"button", value:name})
          .on('mousedown', ()=>{
            $('#processing').css('visibility', '')
          }).on('mouseup', ()=>{
            f();
            $('#processing').css('visibility', 'hidden')
          })
}

function createGenerator(button, type){
  if (!button.hasClass('selected')){
    $('#menu input.selected').removeClass('selected');
    button.addClass('selected');
    let submenu = $('#submenu');
    submenu.empty();
    $("#sheetbar").empty();
    $("#tables").empty();
    let sel;
    switch(type){
      case 'Instructions':
        generateSheet(
            [data.Instructions],
            ["Instructions"],
            ['html'])
        break;
      case 'FileIds':
        generateSheet(
            [data.FileIds],
            ["IDs"],
            ['basic'])
        break;
      case 'HeatMap':
        $('<div/>', {text:"Selected Modification:"}).appendTo(submenu);
        sel = $('<select/>', {name:"mod"}).appendTo(submenu);
        data.HeatMap.map((x, i)=>{
        	$('<option/>', {value:i, text:x.Name}).appendTo(sel);
        })
        longButton("Generate Sheet", ()=>{
          // $("#tables").empty().append(classTableCreate(data.HeatMap[$('#submenu select[name="mod"] option:selected').val()].Data));
          $("#tables").empty()
          generateSheet([data.HeatMap[$('#submenu select[name="mod"] option:selected').val()].Data], [1], ['class'])
        }).appendTo(submenu)
        // $('<input/>', {type:"button", value:"Generate Sheet"}).click().appendTo(submenu);
        // .selectmenu().selectmenu( "menuWidget" ).addClass( "overflow" )
        break;
      case 'ModMapper':
        $('<div/>', {text:"Selected Protein:"}).appendTo(submenu);
        sel = $('<select/>', {name:"protein"}).appendTo(submenu);
        data.ModMapper.map((x, i)=>{
        	$('<option/>', {value:i, text:x.Name}).appendTo(sel);
        })
        longButton("Generate Sheet", ()=>{
          let index = $('#submenu select[name="protein"] option:selected').val(); // index of selected gene
          $("#tables").empty()
          generateSheet(
            [data.ModMapper[index].Summary, ...data.ModMapper[index].Sheets],
            ['Summary', ...[...Array(data.ModMapper[index].Sheets.length).keys()].map(x => ++x)],
            fillArray('basic', data.ModMapper[index].Sheets.length+1)) // all basic sheets
        }).appendTo(submenu)
        // $('<input/>', {type:"button", value:"Generate Sheet"}).click().appendTo(submenu);
        break;
      case 'VennDiagram':
        $('<div/>', {text:"Selected Modification:"}).appendTo(submenu);
        sel = $('<select/>', {name:"mod"}).appendTo(submenu);
        data.VennDiagram.map((x, i)=>{
        	$('<option/>', {value:i, text:x.Name}).appendTo(sel);
        })
        longButton("Generate Sheet", ()=>{
          let index = $('#submenu select[name="mod"] option:selected').val(); // index of selected gene
          $("#tables").empty()

          generateSheet(
            [data.VennDiagram[index].img, ...Object.values(data.VennDiagram[index].raw)],
            ["Image", ...Object.keys(data.VennDiagram[index].raw)],
            ['html', ...fillArray('basic', Object.keys(data.VennDiagram[index].raw).length)]) // all basic sheets
        }).appendTo(submenu)
        break;
      case 'ClusterHeatMap':
	      $('<div/>', {text:"Selected Modification:"}).appendTo(submenu);
        sel = $('<select/>', {name:"mod"}).appendTo(submenu);
        data.ClusterHeatMap.map((x, i)=>{
        	$('<option/>', {value:i, text:x.name}).appendTo(sel);
        })
        longButton("Generate Sheet", ()=>{
          let index = $('#submenu select[name="mod"] option:selected').val(); // index of selected gene
          $("#tables").empty()
          let sheets = data.ClusterHeatMap[index].sheets;

          generateSheet(
            sheets.map(x=>x.data),
            sheets.map(x=>x.name),
            fillArray('html', sheets.length)) // all html sheets
        }).appendTo(submenu)
        break;
      case 'StatTests':
        $('<div/>', {text:"Selected Modification:"}).appendTo(submenu);
        sel = $('<select/>', {name:"mod"}).appendTo(submenu);
        data.StatTests.map((x, i)=>{
          $('<option/>', {value:i, text:x.name}).appendTo(sel);
        })
        longButton("Generate Sheet", ()=>{
          let index = $('#submenu select[name="mod"] option:selected').val(); // index of selected gene
          $("#tables").empty()
          let sheets = data.StatTests[index].data;

          generateSheet(
            Object.values(sheets),
            Object.keys(sheets),
            fillArray('basic', Object.keys(sheets).length)) // all basic sheets
        }).appendTo(submenu)
        break;
      case 'LinearReg':
        $('<div/>', {text:"Selected Modification:"}).appendTo(submenu);
        sel = $('<select/>', {name:"mod"}).appendTo(submenu);
        data.LinearReg.map((x, i)=>{
          $('<option/>', {value:i, text:x.name}).appendTo(sel);
        })
        $('<div/>', {text:"Axis Scale:"}).appendTo(submenu);
        sel = $('<select/>', {name:"axes"}).appendTo(submenu);
        ["norm", "log"].map((x)=>{
          $('<option/>', {value:x, text:x}).appendTo(sel);
        })

        longButton("Generate Sheet", ()=>{
          let index = $('#submenu select[name="mod"] option:selected').val(); // index of selected gene
          let type = $('#submenu select[name="axes"] option:selected').val(); // axes type to select
          $("#tables").empty()
          let sheets = data.LinearReg[index];

          generateSheet(
            [sheets.combined[type], ...sheets.raw.map(x=>x[type])],
            ['combined', ...sheets.raw.map(x=>x.name)],
            fillArray('html', sheets.raw.length+1)) // all basic sheets
        }).appendTo(submenu)
        break;
    }
  }
}

// function to color table
function colorTable(color){
  if (color){
    let nums = [];
    let rlowi = -1, rhighi = -1;
    $('#tables div.container.selected table tr').each(function(index){
      if (index == 0){
        $(this).find('th').each(function(ind){
          let txt = $(this).text();
          if ($.isNumeric(txt)){
            // console.log(txt);
            if(rlowi == -1 && txt === $('input[placeholder="Col Low"]').val()) rlowi = ind;
            if (rhighi == -1 && (txt === $('input[placeholder="Col High"]').val() || $(this).is(':last-child'))) rhighi = ind;
          }
        });
      }else if (index > 1){
        $(this).find('td').each(function(ind){
          if (rlowi <= ind && ind <= rhighi){
            let txt = $(this).text();
            if ($.isNumeric(txt)){
              nums.push(parseInt(txt));
            }
          }
        });
      }
    })
    nums.sort((a, b) => a - b)
    // console.log(nums);
    let pnum = 0; // total number of unique items
    let percentiles = {};
    for(let i = 0; i < nums.length; i++){
      if (!(nums[i] in percentiles)){
        percentiles[nums[i]] = pnum;
        pnum++;
      }
    }
    pnum--;
    // console.log(percentiles);
    // console.log(rlowi, rhighi, max, min);

    $('#tables div.container.selected table tr').each(function(index){
      if (index > 1){
        $(this).find('td').each(function(ind){
          // console.log(rlowi, rhighi);
          if (rlowi <= ind && ind <= rhighi){
            let txt = $(this).text();
            if ($.isNumeric(txt)){
              // this.innerHTML = ind;
              // pickGradient(footerElems.pickrLow._colourValue, footerElems.pickrHigh._colourValue, (parseInt(txt)-min)/max);
              // console.log(Math.floor((footerElems.pickrHigh._colourValue - footerElems.pickrLow._colourValue)*(parseInt(txt)-min)/max));
              $(this).css('background-color', `#${
                pickGradient( footerElems.pickrHigh._colourValue, 
                              footerElems.pickrLow._colourValue, 
                              percentiles[parseInt(txt)]/pnum).toString(16).padStart(6, '0')
              }`);
            }
          }
        });
      }
    })
  }else{
    $('#tables div.container.selected table tr').each(function(index){
      if (index > 1){
        $(this).find('td').each(function(ind){
          $(this).css('background-color', '');
        });
      }
    });
  }
}

var menu = $('#menu');

for(const x of Object.keys(data)){
	$('<input/>', {
	    type:"button",
	    value:x
	    }).click(function(){createGenerator($(this), x)}).appendTo(menu);
}

var footer = $('.box .row.footer');

var footerElems = {};

footerElems.pickrLow = addPickr(footer, "#FFFFFF");
footerElems.pickrLow._root.button.textContent = 'Low Color';
footerElems.pickrHigh = addPickr(footer, "#FF0000");
footerElems.pickrHigh._root.button.textContent = 'High Color';
//Column selectors
$('<input/>', {
  type:"number",
  min:"0",
  placeholder:"Col Low"
  }).change(function(){
    this.value = Math.max(1, Math.floor(this.value));
  }).appendTo(footer);
$('<input/>', {
  type:"number",
  min:"0",
  placeholder:"Col High"
  }).change(function(){
    this.value = Math.max(1, Math.floor(this.value));
  }).appendTo(footer);
//Color tables button
$('<input/>', {
  type:"button",
  value:"Color Table"
  }).click(function(){
    colorTable(true);
  }).appendTo(footer);
//Clear colors button
$('<input/>', {
  type:"button",
  value:"Clear Colors"
  }).click(function(){
    colorTable(false);
  }).appendTo(footer);

//----------------------------------------------------------- MISC UTILITY FUNCTIONS -----------------------------------------------------------

function fillArray(v, n) {
  let a = [];
  for (let i = 0; i < n; i++) {
    a.push(v);
  }
  return a;
}


//colour picker for highlighting
// const initColourStr = "#7289DA";
function addPickr(parent, colourDefault){
  let trigger = document.createElement('button');
  trigger.classList.add('colourPicker');
  trigger.title = 'Set Colour';
  trigger.setAttribute("style", `background-color: ${colourDefault};`);
  parent.append(trigger);
  // $('<button/>', {class:"colorPicker", title:"Set colour", style:`{backgroundColor:${initColourStr}}`})
  let ret = new Pickr({
          el: trigger,
          theme: 'monolith',
          useAsButton: true,
          default: colourDefault,
          position: 'top-middle',
          components: {
            preview: true,
            hue: true,
            interaction: {
              input: true,
              clear: true
            },
          }
        }).on('change', colour => {
          trigger.style.backgroundColor = colour.toRGBA()
        }).on('hide', instance => {
          instance._colourValue = parseInt(instance.getColor().toHEXA().toString().slice(1), 16)
          // if (window.onChange) window.onChange()
        }).on('clear', instance => {
          trigger.style.backgroundColor = instance.options.default;
          instance._colourValue = parseInt(instance.options.default.slice(1), 16)
        })
  ret._colourValue = parseInt(colourDefault.slice(1), 16);
  return ret;
}