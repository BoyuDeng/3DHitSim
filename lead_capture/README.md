# Lead Capture Software

A lightweight lead capture web app you can run locally with no external dependencies.

## Features

- Web form for capturing leads (name, email, phone, company, source, notes, consent)
- SQLite persistence
- Recent leads table on the homepage
- JSON API endpoint (`/api/leads`)
- CSV export endpoint (`/export.csv`)

## Run

```bash
python3 lead_capture/server.py --host 0.0.0.0 --port 8080 --db lead_capture/leads.db
```

Then open: <http://localhost:8080>

## Endpoints

- `GET /` - form + recent leads
- `POST /submit` - save a lead
- `GET /api/leads` - list leads as JSON
- `GET /export.csv` - download leads as CSV
