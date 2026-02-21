import json
import tempfile
import threading
import time
import unittest
from pathlib import Path
from urllib.parse import urlencode
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from lead_capture.server import make_handler, LeadStore
from http.server import ThreadingHTTPServer


class LeadCaptureTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        db = Path(self.tmp.name) / "test.db"
        store = LeadStore(db)
        self.server = ThreadingHTTPServer(("127.0.0.1", 0), make_handler(store))
        self.port = self.server.server_address[1]
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        time.sleep(0.05)

    def tearDown(self):
        self.server.shutdown()
        self.server.server_close()
        self.thread.join(timeout=2)
        self.tmp.cleanup()

    def url(self, path):
        return f"http://127.0.0.1:{self.port}{path}"

    def test_home_page(self):
        with urlopen(self.url("/")) as resp:
            body = resp.read().decode("utf-8")
        self.assertIn("Lead Capture", body)

    def test_submit_and_api(self):
        data = urlencode({"name": "Alice", "email": "alice@example.com", "consent": "on"}).encode()
        req = Request(self.url("/submit"), data=data, method="POST")
        with urlopen(req) as resp:
            self.assertEqual(resp.status, 200)

        with urlopen(self.url("/api/leads")) as resp:
            payload = json.loads(resp.read().decode("utf-8"))
        self.assertEqual(payload["count"], 1)
        self.assertEqual(payload["leads"][0]["name"], "Alice")

    def test_validation(self):
        data = urlencode({"name": "", "email": "bad"}).encode()
        req = Request(self.url("/submit"), data=data, method="POST")
        with self.assertRaises(HTTPError) as ctx:
            urlopen(req)
        self.assertEqual(ctx.exception.code, 400)
        body = ctx.exception.read().decode("utf-8")
        self.assertIn("Name is required.", body)
        self.assertIn("Email format is invalid.", body)


if __name__ == "__main__":
    unittest.main()
