import csv, json, pathlib, subprocess, sys, tempfile, unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
class PipelineTest(unittest.TestCase):
    def test_example(self):
        with tempfile.TemporaryDirectory() as tmp:
            subprocess.run([sys.executable,str(ROOT/'full_pipeline.py'),'--input',str(ROOT/'data'),'--output',tmp],check=True,capture_output=True)
            with open(pathlib.Path(tmp)/'resolved_alleles.csv') as f: alleles=list(csv.DictReader(f))
            with open(pathlib.Path(tmp)/'family_diagnoses.csv') as f: diag=list(csv.DictReader(f))
            qc=json.loads((pathlib.Path(tmp)/'qc_report.json').read_text())
            self.assertEqual(len(alleles),190)
            self.assertEqual(len(diag),32)
            self.assertEqual(qc['cnv_merges'],3)
            self.assertTrue(set(d['status'] for d in diag)<= {'diagnosed','negative','unresolved'})
if __name__=='__main__':unittest.main()
