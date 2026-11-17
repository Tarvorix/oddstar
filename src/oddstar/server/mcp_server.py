"""Model Context Protocol (MCP) tool-server wrapper over the oddstar CLI (Section 10).

The CLI remains the single source of truth; every MCP tool constructs the exact
CLI invocation and returns the structured JSON contract.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

from oddstar.cli import dispatch
from oddstar.config import repo_root

MCP_TOOLS: list[dict[str, Any]] = [
    {
        "name": "oddstar_doctor",
        "description": "Check Python environment, lock file, cache cap, and archive definitions.",
        "inputSchema": {"type": "object", "properties": {"release": {"type": "string", "enum": ["dr3", "dr4"]}}},
    },
    {
        "name": "oddstar_sources_list",
        "description": "List all 23 registered public archives, releases, and rate limits.",
        "inputSchema": {"type": "object", "properties": {}},
    },
    {
        "name": "oddstar_fetch",
        "description": "Fetch raw archive data into content-addressed cache and write manifest.",
        "inputSchema": {
            "type": "object",
            "required": ["archive"],
            "properties": {
                "archive": {"type": "string"},
                "ids": {"type": "string"},
                "cone": {"type": "array", "items": {"type": "number"}, "minItems": 3, "maxItems": 3},
                "adql": {"type": "string"},
                "dry_run": {"type": "boolean"},
            },
        },
    },
    {
        "name": "oddstar_sample_build",
        "description": "Build parent star sample for a channel with cut accounting.",
        "inputSchema": {
            "type": "object",
            "required": ["channel"],
            "properties": {"channel": {"type": "string"}, "limit": {"type": "integer"}},
        },
    },
    {
        "name": "oddstar_validate",
        "description": "Run 4-part validation harness and issue a frozen certificate.",
        "inputSchema": {
            "type": "object",
            "required": ["channel"],
            "properties": {
                "channel": {"type": "string"},
                "split": {"type": "string", "enum": ["dev", "holdout", "all"]},
                "release": {"type": "string", "enum": ["dr3", "dr4"]},
                "limit": {"type": "integer"},
            },
        },
    },
    {
        "name": "oddstar_search",
        "description": "Run certified channel search (blocked without valid certificate).",
        "inputSchema": {
            "type": "object",
            "required": ["channel"],
            "properties": {
                "channel": {"type": "string"},
                "run_id": {"type": "string"},
                "limit": {"type": "integer"},
                "seed": {"type": "integer"},
                "extra_ids": {"type": "string"},
            },
        },
    },
    {
        "name": "oddstar_vet",
        "description": "Run tri-state artifact and astrophysical impostor checks on flagged targets.",
        "inputSchema": {"type": "object", "required": ["run_id"], "properties": {"run_id": {"type": "string"}}},
    },
    {
        "name": "oddstar_fuse",
        "description": "Join channel results on source_id, calibrate FAP, assign tiers, and rank.",
        "inputSchema": {
            "type": "object",
            "required": ["run_ids"],
            "properties": {"run_ids": {"type": "array", "items": {"type": "string"}}, "out_run_id": {"type": "string"}},
        },
    },
    {
        "name": "oddstar_dossier",
        "description": "Generate dossier.json, dossier.md, and evidence plots for a source_id.",
        "inputSchema": {"type": "object", "required": ["source_id"], "properties": {"source_id": {"type": "string"}}},
    },
    {
        "name": "oddstar_inspect",
        "description": "One-star full inspection mode across channels, vetting, and dossier.",
        "inputSchema": {
            "type": "object",
            "required": ["source_id"],
            "properties": {"source_id": {"type": "string"}, "release": {"type": "string", "enum": ["dr3", "dr4"]}},
        },
    },
    {
        "name": "oddstar_review_packet",
        "description": "Build blinded reviewer packet mixed with controls and known objects.",
        "inputSchema": {
            "type": "object",
            "required": ["ids"],
            "properties": {"ids": {"type": "string"}, "round": {"type": "integer"}, "claim": {"type": "string"}},
        },
    },
    {
        "name": "oddstar_review_ingest",
        "description": "Validate and ingest reviewer JSON record.",
        "inputSchema": {"type": "object", "required": ["file"], "properties": {"file": {"type": "string"}}},
    },
    {
        "name": "oddstar_lessons_add",
        "description": "Record a lesson closing an open review finding with a linked test.",
        "inputSchema": {
            "type": "object",
            "required": [
                "lesson_id",
                "finding_id",
                "source_id",
                "title",
                "what_happened",
                "generic_lesson",
                "automated_check",
                "linked_test",
            ],
            "properties": {
                "lesson_id": {"type": "string"},
                "finding_id": {"type": "string"},
                "source_id": {"type": "string"},
                "title": {"type": "string"},
                "what_happened": {"type": "string"},
                "generic_lesson": {"type": "string"},
                "automated_check": {"type": "string"},
                "linked_test": {"type": "string"},
            },
        },
    },
    {
        "name": "oddstar_null_report",
        "description": "Emit sensitivity completeness report and upper limits for a run.",
        "inputSchema": {"type": "object", "required": ["run_id"], "properties": {"run_id": {"type": "string"}, "channel": {"type": "string"}}},
    },
    {
        "name": "oddstar_status",
        "description": "Report certified channels, active runs, and open review findings.",
        "inputSchema": {"type": "object", "properties": {}},
    },
    {
        "name": "oddstar_show",
        "description": "Inspect JSON/JSONL/text artifacts in bounded pages.",
        "inputSchema": {
            "type": "object",
            "required": ["path"],
            "properties": {"path": {"type": "string"}, "rows": {"type": "integer"}, "offset": {"type": "integer"}},
        },
    },
    {
        "name": "oddstar_schema",
        "description": "Print or export JSON Schemas.",
        "inputSchema": {
            "type": "object",
            "properties": {"type": {"type": "string"}, "export": {"type": "boolean"}},
        },
    },
]


def _tool_args_to_argv(tool_name: str, arguments: dict[str, Any]) -> list[str]:
    argv: list[str] = []
    if tool_name == "oddstar_doctor":
        argv = ["doctor"]
    elif tool_name == "oddstar_sources_list":
        argv = ["sources", "list"]
    elif tool_name == "oddstar_fetch":
        argv = ["fetch", str(arguments["archive"])]
        if arguments.get("ids"):
            argv.extend(["--ids", str(arguments["ids"])])
        if arguments.get("cone"):
            argv.extend(["--cone", *[str(x) for x in arguments["cone"]]])
        if arguments.get("adql"):
            argv.extend(["--adql", str(arguments["adql"])])
    elif tool_name == "oddstar_sample_build":
        argv = ["sample", "build", str(arguments["channel"])]
    elif tool_name == "oddstar_validate":
        argv = ["validate", str(arguments["channel"])]
        if arguments.get("split"):
            argv.extend(["--split", str(arguments["split"])])
    elif tool_name == "oddstar_search":
        argv = ["search", str(arguments["channel"])]
        if arguments.get("run_id"):
            argv.extend(["--run-id", str(arguments["run_id"])])
        if arguments.get("seed") is not None:
            argv.extend(["--seed", str(arguments["seed"])])
        if arguments.get("extra_ids"):
            argv.extend(["--extra-ids", str(arguments["extra_ids"])])
    elif tool_name == "oddstar_vet":
        argv = ["vet", str(arguments["run_id"])]
    elif tool_name == "oddstar_fuse":
        argv = ["fuse", *[str(r) for r in arguments["run_ids"]]]
        if arguments.get("out_run_id"):
            argv.extend(["--out-run-id", str(arguments["out_run_id"])])
    elif tool_name == "oddstar_dossier":
        argv = ["dossier", str(arguments["source_id"])]
    elif tool_name == "oddstar_inspect":
        argv = ["inspect", str(arguments["source_id"])]
    elif tool_name == "oddstar_review_packet":
        argv = ["review", "packet", str(arguments["ids"])]
        if arguments.get("round") is not None:
            argv.extend(["--round", str(arguments["round"])])
        if arguments.get("claim"):
            argv.extend(["--claim", str(arguments["claim"])])
    elif tool_name == "oddstar_review_ingest":
        argv = ["review", "ingest", str(arguments["file"])]
    elif tool_name == "oddstar_lessons_add":
        argv = [
            "lessons",
            "add",
            "--lesson-id",
            str(arguments["lesson_id"]),
            "--finding-id",
            str(arguments["finding_id"]),
            "--source-id",
            str(arguments["source_id"]),
            "--title",
            str(arguments["title"]),
            "--what-happened",
            str(arguments["what_happened"]),
            "--generic-lesson",
            str(arguments["generic_lesson"]),
            "--automated-check",
            str(arguments["automated_check"]),
            "--linked-test",
            str(arguments["linked_test"]),
        ]
    elif tool_name == "oddstar_null_report":
        argv = ["null-report", str(arguments["run_id"])]
        if arguments.get("channel"):
            argv.extend(["--channel", str(arguments["channel"])])
    elif tool_name == "oddstar_status":
        argv = ["status"]
    elif tool_name == "oddstar_show":
        argv = ["show", str(arguments["path"])]
        if arguments.get("rows") is not None:
            argv.extend(["--rows", str(arguments["rows"])])
        if arguments.get("offset") is not None:
            argv.extend(["--offset", str(arguments["offset"])])
    elif tool_name == "oddstar_schema":
        argv = ["schema"]
        if arguments.get("type"):
            argv.append(str(arguments["type"]))
        if arguments.get("export"):
            argv.append("--export")
    else:
        raise ValueError(f"Unknown MCP tool: {tool_name}")

    if arguments.get("limit") is not None:
        argv.extend(["--limit", str(arguments["limit"])])
    if arguments.get("release"):
        argv.extend(["--release", str(arguments["release"])])
    if arguments.get("dry_run"):
        argv.append("--dry-run")
    argv.append("--json")
    return argv


def call_tool(tool_name: str, arguments: dict[str, Any], *, root: Path | None = None) -> tuple[int, dict[str, Any]]:
    """Invoke an oddstar MCP tool using the CLI contract."""
    argv = _tool_args_to_argv(tool_name, arguments)
    return dispatch(argv, root=root or repo_root())


def handle_jsonrpc_request(req: dict[str, Any], *, root: Path | None = None) -> dict[str, Any]:
    """Handle one MCP JSON-RPC 2.0 request."""
    req_id = req.get("id")
    method = req.get("method")
    params = req.get("params") or {}

    if method == "initialize":
        return {
            "jsonrpc": "2.0",
            "id": req_id,
            "result": {
                "protocolVersion": "2024-11-05",
                "serverInfo": {"name": "oddstar-mcp", "version": "0.1.0"},
                "capabilities": {"tools": {}},
            },
        }
    if method == "tools/list":
        return {"jsonrpc": "2.0", "id": req_id, "result": {"tools": MCP_TOOLS}}
    if method == "tools/call":
        name = params.get("name", "")
        args = params.get("arguments") or {}
        code, payload = call_tool(name, args, root=root)
        return {
            "jsonrpc": "2.0",
            "id": req_id,
            "result": {
                "content": [{"type": "text", "text": json.dumps(payload, indent=2)}],
                "isError": code != 0,
                "exitCode": code,
            },
        }
    return {
        "jsonrpc": "2.0",
        "id": req_id,
        "error": {"code": -32601, "message": f"Method not found: {method}"},
    }


def serve_stdio() -> int:
    """Run line-delimited JSON-RPC 2.0 MCP stdio server."""
    for line in sys.stdin:
        if not line.strip():
            continue
        req = json.loads(line)
        resp = handle_jsonrpc_request(req)
        sys.stdout.write(json.dumps(resp) + "\n")
        sys.stdout.flush()
    return 0


if __name__ == "__main__":
    raise SystemExit(serve_stdio())
