importScripts('call_gist.js');
onmessage = (e) => {
    console.log("zoombox worker "+self.name);
    var result = (Module.do_forbid_w(e.data[0],e.data[1],e.data[2],e.data[3],e.data[4],e.data[5],e.data[6],e.data[7],e.data[8],e.data[9], e.data[10],e.data[11],e.data[12],e.data[13],e.data[14]));
    result = JSON.parse(result);
    postMessage({'result':result,'id':e.data[15]});
    console.log("finished job worker "+self.name);
}
