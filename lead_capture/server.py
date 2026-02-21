#!/usr/bin/env python3
"""Simple lead capture web app using only Python stdlib."""

from __future__ import annotations

import argparse
import csv
import html
import io
import json
import re
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Dict, Iterable, List, Tuple
from urllib.parse import parse_qs, urlparse

EMAIL_RE = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")


@dataclass
class Lead:
    name: str
    email: str
    phone: str = ""
    company: str = ""
    source: str = ""
    notes: str = ""
    consent: bool = False


class LeadStore:
    def __init__(self, db_path: Path) -> None:
        self.db_path = db_path
        self._init_db()

    def _connect(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        return conn

    def _init_db(self) -> None:
        with self._connect() as conn:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS leads (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    name TEXT NOT NULL,
                    email TEXT NOT NULL,
                    phone TEXT,
                    company TEXT,
                    source TEXT,
                    notes TEXT,
                    consent INTEGER NOT NULL DEFAULT 0,
                    created_at TEXT NOT NULL
                )
                """
            )

    def add_lead(self, lead: Lead) -> int:
        with self._connect() as conn:
            cur = conn.execute(
                """
                INSERT INTO leads (
                    name, email, phone, company, source, notes, consent, created_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    lead.name,
                    lead.email,
                    lead.phone,
                    lead.company,
                    lead.source,
                    lead.notes,
                    int(lead.consent),
                    datetime.now(timezone.utc).isoformat(),
                ),
            )
            return int(cur.lastrowid)

    def list_leads(self, limit: int = 50) -> List[Dict[str, object]]:
        with self._connect() as conn:
            rows = conn.execute(
                """
                SELECT id, name, email, phone, company, source, notes, consent, created_at
                FROM leads
                ORDER BY id DESC
                LIMIT ?
                """,
                (limit,),
            ).fetchall()
        return [dict(r) for r in rows]


class LeadCaptureHandler(BaseHTTPRequestHandler):
    store: LeadStore

    def do_GET(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)
        if parsed.path == "/":
            self._serve_home()
            return
        if parsed.path == "/api/leads":
            self._serve_leads_json()
            return
        if parsed.path == "/export.csv":
            self._serve_csv()
            return
        self.send_error(HTTPStatus.NOT_FOUND, "Not Found")

    def do_POST(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)
        if parsed.path != "/submit":
            self.send_error(HTTPStatus.NOT_FOUND, "Not Found")
            return

        length = int(self.headers.get("Content-Length", "0"))
        raw = self.rfile.read(length).decode("utf-8")
        form = {k: v[0] if v else "" for k, v in parse_qs(raw).items()}

        errors = self._validate(form)
        if errors:
            self._serve_home(errors=errors, form=form, status=HTTPStatus.BAD_REQUEST)
            return

        lead = Lead(
            name=form.get("name", "").strip(),
            email=form.get("email", "").strip(),
            phone=form.get("phone", "").strip(),
            company=form.get("company", "").strip(),
            source=form.get("source", "").strip(),
            notes=form.get("notes", "").strip(),
            consent=form.get("consent") == "on",
        )
        self.store.add_lead(lead)
        self._redirect("/?saved=1")

    def _validate(self, form: Dict[str, str]) -> List[str]:
        errors: List[str] = []
        if not form.get("name", "").strip():
            errors.append("Name is required.")
        email = form.get("email", "").strip()
        if not email:
            errors.append("Email is required.")
        elif not EMAIL_RE.match(email):
            errors.append("Email format is invalid.")
        return errors

    def _serve_home(
        self,
        errors: Iterable[str] = (),
        form: Dict[str, str] | None = None,
        status: HTTPStatus = HTTPStatus.OK,
    ) -> None:
        form = form or {}
        saved = "saved=1" in urlparse(self.path).query
        leads = self.store.list_leads(20)

        error_html = "".join(f"<li>{html.escape(e)}</li>" for e in errors)
        flash = (
            "<p class='success'>Lead saved successfully.</p>" if saved else ""
        ) + (f"<ul class='errors'>{error_html}</ul>" if error_html else "")

        def val(name: str) -> str:
            return html.escape(form.get(name, ""))

        rows = "".join(
            "<tr>"
            f"<td>{lead['id']}</td>"
            f"<td>{html.escape(str(lead['name']))}</td>"
            f"<td>{html.escape(str(lead['email']))}</td>"
            f"<td>{html.escape(str(lead['phone'] or ''))}</td>"
            f"<td>{html.escape(str(lead['company'] or ''))}</td>"
            f"<td>{html.escape(str(lead['source'] or ''))}</td>"
            f"<td>{'yes' if lead['consent'] else 'no'}</td>"
            f"<td>{html.escape(str(lead['created_at']))}</td>"
            "</tr>"
            for lead in leads
        )

        page = f"""<!doctype html>
<html>
<head>
  <meta charset='utf-8' />
  <title>Lead Capture</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 2rem; max-width: 980px; }}
    form {{ display: grid; grid-template-columns: repeat(2, 1fr); gap: 0.8rem; margin-bottom: 1rem; }}
    .full {{ grid-column: 1 / -1; }}
    input, textarea {{ width: 100%; padding: 0.4rem; }}
    button {{ padding: 0.6rem 1rem; }}
    table {{ border-collapse: collapse; width: 100%; }}
    th, td {{ border: 1px solid #ddd; padding: 0.45rem; text-align: left; }}
    .errors {{ color: #a00; }}
    .success {{ color: #067d17; }}
    .top-actions {{ margin: 0.8rem 0; }}
  </style>
</head>
<body>
  <h1>Lead Capture</h1>
  {flash}
  <form method='POST' action='/submit'>
    <label>Name*<br><input type='text' name='name' value='{val('name')}' required></label>
    <label>Email*<br><input type='email' name='email' value='{val('email')}' required></label>
    <label>Phone<br><input type='text' name='phone' value='{val('phone')}'></label>
    <label>Company<br><input type='text' name='company' value='{val('company')}'></label>
    <label class='full'>Source<br><input type='text' name='source' value='{val('source')}' placeholder='Landing page, ad campaign, referral, etc.'></label>
    <label class='full'>Notes<br><textarea name='notes' rows='4'>{val('notes')}</textarea></label>
    <label class='full'><input type='checkbox' name='consent' {'checked' if form.get('consent') == 'on' else ''}> I consent to follow-up communications</label>
    <div class='full'><button type='submit'>Save Lead</button></div>
  </form>

  <div class='top-actions'>
    <a href='/api/leads'>View JSON API</a> |
    <a href='/export.csv'>Export CSV</a>
  </div>

  <h2>Recent Leads</h2>
  <table>
    <thead><tr><th>ID</th><th>Name</th><th>Email</th><th>Phone</th><th>Company</th><th>Source</th><th>Consent</th><th>Created (UTC)</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
</body>
</html>"""
        body = page.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _serve_leads_json(self) -> None:
        leads = self.store.list_leads(500)
        payload = json.dumps({"count": len(leads), "leads": leads}, indent=2).encode("utf-8")
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def _serve_csv(self) -> None:
        leads = self.store.list_leads(5000)
        output = io.StringIO()
        writer = csv.DictWriter(
            output,
            fieldnames=["id", "name", "email", "phone", "company", "source", "notes", "consent", "created_at"],
        )
        writer.writeheader()
        for lead in leads:
            writer.writerow(lead)
        data = output.getvalue().encode("utf-8")

        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "text/csv; charset=utf-8")
        self.send_header("Content-Disposition", 'attachment; filename="leads.csv"')
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def _redirect(self, location: str) -> None:
        self.send_response(HTTPStatus.SEE_OTHER)
        self.send_header("Location", location)
        self.end_headers()


def make_handler(store: LeadStore):
    class _Handler(LeadCaptureHandler):
        pass

    _Handler.store = store
    return _Handler


def run(host: str, port: int, db_path: Path) -> None:
    store = LeadStore(db_path)
    server = ThreadingHTTPServer((host, port), make_handler(store))
    print(f"Lead capture server running at http://{host}:{port} (db: {db_path})")
    server.serve_forever()


def main() -> None:
    parser = argparse.ArgumentParser(description="Run lead capture software")
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--port", type=int, default=8080)
    parser.add_argument("--db", type=Path, default=Path("leads.db"))
    args = parser.parse_args()
    run(args.host, args.port, args.db)


if __name__ == "__main__":
    main()
