"""MCP tool server wrapping the oddstar CLI."""

from oddstar.server.mcp_server import MCP_TOOLS, call_tool, handle_jsonrpc_request

__all__ = ["MCP_TOOLS", "call_tool", "handle_jsonrpc_request"]
