import csv, json, pathlib, subprocess, sys, tempfile, unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
class PipelineTest(unittest.TestCase):
    def test_small_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            subprocess.run([sys.executable,str(ROOT/'full_pipeline.py'),'--input',str(ROOT/'data'),'--output',tmp,'--limit','3'],check=True,capture_output=True)
            with open(pathlib.Path(tmp)/'exposure.csv') as f: exposure=list(csv.DictReader(f))
            with open(pathlib.Path(tmp)/'regimen_scores.csv') as f: score=list(csv.DictReader(f))
            qc=json.loads((pathlib.Path(tmp)/'qc_report.json').read_text())
            self.assertEqual(len(exposure),3)
            self.assertEqual(len(score),14)
            self.assertEqual(qc['courses_processed'],3)
            self.assertEqual(sorted(set(int(r['rank']) for r in score)),list(range(1,8)))
if __name__=='__main__':unittest.main()
