load ("extras/mongrel2.jl")

t = m2_run_server("6DFF1523-C091-49B8-B635-598640E864B3", "tcp://127.0.0.1:9997", "tcp://127.0.0.1:9996")
while true
	(conn, req) = consume (t) 	
	println("$(req.json_data) $(req.connection_id) $(req.sender_id)"); flush(stdout_stream);

	if m2_is_disconnected(req); println("Disconnected $(req.connection_id)"); continue; end

	m2_reply(conn, req, "{\"queue\":[[1]]}")
end
