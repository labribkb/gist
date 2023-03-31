//importScripts('call_forbid.js');
const pool = new Map();
for(let i=0;i<2;i++){
    const name = "worker_from_pool_"+i;
    const worker = new Worker("workerZoomBox.js",{"name":name});
    console.log(worker);
    pool.set(worker,0);
    worker.onmessage =((e) => {
        if(e.data['result'] !==undefined){
            pool.set(worker,pool.get(worker)-1);
            console.log("a worker has finished "+worker);
        }
        postMessage(e.data);
    });
}
onmessage = (e) => {
    if(typeof(e.data)=="string" && e.data=="terminate"){
        for (const [worker, nbJobs] of pool) {
            worker.terminate();
        }
    }else{
        let found = false;
        let minJobs= pool[Symbol.iterator]().next().value;
        for (const [worker, nbJobs] of pool) {
            if(nbJobs==0){
                pool.set(worker,1);
                console.log("choosing worker "+worker);
                worker.postMessage(e.data);
                found=true;
                break;
            }else{
                if(minJobs>nbJobs){
                    minJobs=nbJobs;
                }
            }
        }
        if(!found){
            for (const [worker, nbJobs] of pool) {
                if(nbJobs==minJobs){
                    pool.set(worker,nbJobs+1);
                    console.log("choosing worker "+worker);
                    worker.postMessage(e.data);
                    break;
                }
            }
        }
    }
}
