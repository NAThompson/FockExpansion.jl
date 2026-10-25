"""Run one scientific script in-process with periodic RSS and CPU guards."""
import sys
sys.dont_write_bytecode=True
import os,time,json,resource,signal,runpy,contextlib,traceback,fcntl
from pathlib import Path
out=Path(__file__).resolve().parent
name=sys.argv[1]
assert name.startswith('x24_') and '/' not in name and name.endswith('.py')
stem=name[:-3]; start=time.time(); rc=0;reason=None;peak=0
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
os.environ['TMPDIR']=str(out)
resource.setrlimit(resource.RLIMIT_CPU,(600,600))
def guard(signum,frame):
    global peak
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    val=rss if sys.platform=='darwin' else rss*1024
    peak=max(peak,val)
    if val>1_000_000_000:raise MemoryError('x24 1 GB RSS stop threshold')
    if time.time()-start>650:raise TimeoutError('x24 650 s wall-clock stop')
with (out/'x24_science.lock').open('a') as lock:
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.signal(signal.SIGALRM,guard);signal.setitimer(signal.ITIMER_REAL,.1,.1)
    with (out/(stem+'.out')).open('a') as stdout,(out/(stem+'.err')).open('a') as stderr:
        with contextlib.redirect_stdout(stdout),contextlib.redirect_stderr(stderr):
            try:runpy.run_path(str(out/name),run_name='__main__')
            except BaseException as exc:
                traceback.print_exc();rc=1;reason=type(exc).__name__+': '+str(exc)
    signal.setitimer(signal.ITIMER_REAL,0)
    guard(None,None)
info={'script':name,'exit_code':rc,'reason':reason,'peak_rss_bytes':peak,'elapsed_seconds':time.time()-start,
      'guard':'single scientific process; SIGALRM checks ru_maxrss every 0.1 s and stops at 1 GB; CPU limit 600 s; wall stop 650 s; exclusive scientific lock',
      'limits_note':'Darwin rejects lowering RLIMIT_AS/DATA/RSS; sandbox forbids ps. RSS guard runs within the scientific process and does not spawn a monitoring job.'}
(out/(stem+'_resource.json')).write_text(json.dumps(info,indent=2)+'\n')
print(json.dumps(info,indent=2));sys.exit(rc)
