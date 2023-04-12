function load_file(inputData,inputLabels){
	var fileReader = new FileReader();
	fileReader.readAsText(inputData);
	fileReader.onload= function(){
        if(inputLabels.length>0){
            var fileReaderLabels = new FileReader();
            fileReaderLabels.readAsText(inputLabels);
            fileReaderLabels.onload= function(){
                Module.load(fileReader.result,fileReaderLabels.result);
            };
            fileReaderLabels.onerror= function(){
                alert(fileReader.error);
            }
        }else{
                Module.load(fileReader.result,"");
        }
	}
	fileReader.onerror= function(){
		alert(fileReader.error);
	}
}
function startload(){
    const inputData =document.getElementById("dataFile");
    const inputLabels =document.getElementById("labelFile");
    if(inputData.files.length==1&&inputLabels.files.length==1)
        load_file(inputData.files[0],inputLabels.files[0]);
}  

function startloadfromlist(){
    const inputData =document.getElementById("datasetchoice");

    Module.loadFromDataList(parseInt(inputData.value));
    $("#datasetdescription dl").addClass("hiddendesc");
    $("#datasetdescription dl:nth-child("+(1+parseInt(inputData.value))+")").removeClass("hiddendesc");
}
$(function(){
    $( "#accordion" ).accordion();
    $(document).uitooltip({
        items:"#instructionrect",
        content:function(){
            return "<img class='instructiontooltip' src='instructionrect.gif'>";
        }
    });
    $("#tolerance").slider({
        range:"min",
        min:0,
        max:32,
        value:3,
        slide:function(event,ui){
            $("#tolerancetext").val(ui.value);
            Module.changeTolerance($("#tolerance").slider("value"),false);
        }
    });
    $("#tolerancetext").val($("#tolerance").slider("value").toString());
    $("#updatetolerance").button();
    $("#updatetolerance").on("click",function(e){
        Module.changeTolerance($("#tolerance").slider("value"),true);
    });
    $('#download-choice').dialog({
        resizable:false,
        autoOpen:false,
        buttons:[
            {
                text:"As PNG Image",
                click: function() {
                    Module.download_render();
                    $(this).dialog("close");
                }
            },{
                text:"As Text file",
                click: function() {
                    Module.download_layout();
                    $(this).dialog("close");
                }
            }
        ],
        modal:true
    });
    $('#download').on('click',function(){
        $("#download-choice").dialog( "open" );
    });
});
function onzoombox(n,dataurl,w,h){
    $("<h3>Zoom on "+n+" elements</h3><div><img src=\""+dataurl+"\" width="+w+" height="+h+" ></div>").appendTo('#accordion');
    $( "#accordion" ).accordion("refresh");
    $( "#accordion" ).accordion({active:-1});
}
function onunzoombox(){
    $( "#accordion>:last-of-type" ).remove();
    $( "#accordion" ).accordion("refresh");
    $( "#accordion" ).accordion({active:-1});
}
function adjust_res_message(R){
    $('#res_message>#res').text(`${R} x ${R}`);
    if(R==1024){
        $('#res_message').css("visibility","hidden");
    }else{
        $('#res_message').css("visibility","unset");
    }
}
function hide_best_message(){
    $('#best_params').css('visibility','hidden');
}
function adjust_best_message(R,best_tol,best_diam){
    var d = parseFloat(best_diam.toFixed(2));
    var t = parseFloat(best_tol.toFixed(2));
    $('#best_params').css('visibility','unset');
    $('#best_params #res').text(`${R} x ${R}`);
    $('#best_params #best_diam').text(`${d} pixels in visual space`);
    $('#best_params #best_tol').text(`${t} pixels in visual space`);
}
function download_canvas_image(){
    let l = document.createElement('a');
    l.href = Module.canvas.toDataURL();
    l.download = 'gist.png';
    l.click();
}

function download_layout(layout){
    let l = document.createElement('a');
    l.href = `data:text/plain,${layout}`;
    l.download = 'gist.txt';
    l.click();
}

