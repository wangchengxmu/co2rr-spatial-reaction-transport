"""Portable runner contract tests; no nonlinear simulations launched."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace
import run_screen as r

class ScreenTests(unittest.TestCase):
    def test_grid(self):
        rows=r.make_cases();self.assertEqual(len(rows),30)
        self.assertEqual(len({x['id'] for x in rows}),30)
        self.assertEqual({x['delta'] for x in rows},{10,30,50})
        self.assertEqual({x['chi'] for x in rows},set(r.CHIS))
        self.assertTrue(all((x['nc'],x['ny'],x['nx'],x['outer'])==(48,96,24,500) for x in rows))
        all_rows=r.make_cases(True);self.assertEqual(len(all_rows),36)
        self.assertTrue(all(x['chi']==.05 for x in all_rows[30:]))

    def test_model_and_process_must_both_pass(self):
        for model_status,code,expected in [('passed',0,'passed'),('not_accepted',0,'not_accepted'),('passed',1,'not_accepted')]:
            with self.subTest(status=model_status,code=code),tempfile.TemporaryDirectory() as tmp:
                out=Path(tmp);case=r.make_cases()[0]
                def worker(command,**kwargs):
                    target=Path(command[command.index('--out')+1])
                    (target/'result.json').write_text(json.dumps({'status':model_status}),encoding='utf-8')
                    return SimpleNamespace(returncode=code)
                with patch.object(r.subprocess,'run',side_effect=worker):
                    actual=r._run_case(case,Path(__file__).resolve().parent,out)
                self.assertEqual(actual['status'],expected)

    def test_refuses_nonempty_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);(out/'existing.txt').write_text('preserve')
            with self.assertRaises(FileExistsError):
                r.run(SimpleNamespace(out=str(out),verification=False,limit=1,workers=1))
            self.assertEqual((out/'existing.txt').read_text(),'preserve')

if __name__=='__main__':unittest.main()
