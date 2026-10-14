"""Run the notebook cell by cell, saving progress and a readable log.
Usage: .venv/bin/python run_notebook.py [first_cell] [notebook.ipynb]"""
import sys, time, nbformat
from nbclient import NotebookClient
nbfile = sys.argv[2] if len(sys.argv) > 2 else 'Capra_hircus_analysis.ipynb'
start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
out_file = nbfile.replace('.ipynb', '_executed.ipynb')
nb = nbformat.read(nbfile, as_version=4)
nb.metadata['kernelspec'] = {'name': 'capra', 'display_name': 'Capra hircus', 'language': 'python'}
client = NotebookClient(nb, timeout=None, kernel_name='capra', resources={'metadata': {'path': '.'}})
log = open('run.log', 'a')
def say(m): log.write(m + '\n'); log.flush()
with client.setup_kernel():
    for i, cell in enumerate(nb.cells):
        if cell.cell_type != 'code' or i < start: continue
        t0 = time.time(); say(f'--- cell {i}: {cell.source.splitlines()[0][:70]}')
        try:
            client.execute_cell(cell, i)
        except Exception as e:
            say(f'ERROR in cell {i}: {type(e).__name__}')
            for o in cell.get('outputs', []):
                say(o.get('text') or '\n'.join(o.get('traceback', [])) if o.output_type in ('stream', 'error') else '')
            nbformat.write(nb, out_file); sys.exit(1)
        for o in cell.get('outputs', []):
            if o.output_type == 'stream': say(o.text[-1500:])
            elif o.output_type in ('execute_result', 'display_data') and 'text/plain' in o.data: say(o.data['text/plain'][-1500:])
        say(f'    ok ({time.time()-t0:.0f} s)')
        nbformat.write(nb, out_file)
say('DONE: all cells executed')
