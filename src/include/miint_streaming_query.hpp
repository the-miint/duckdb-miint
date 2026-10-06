#pragma once

#include "duckdb/main/connection.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/main/query_result_stream.hpp"

namespace duckdb {

//! A streaming query result with the interface miint relied on from v1.5's Connection::SendQuery.
//! DuckDB v2.0 replaced SendQuery with Submit + QueryResultStream. A QueryResultStream cannot be opened on a result
//! that already failed (bind or plan error) -- its constructor throws -- so the error is kept here and reported via
//! HasError()/GetError(), letting callers keep their own messages. Execution errors surface as in v1.5: Fetch()
//! returns nullptr and HasError() becomes true. Fetched chunks are flat (ChunkFormat copies them for buffering).
class StreamingQuery {
public:
	explicit StreamingQuery(unique_ptr<QueryResult> submitted) {
		if (submitted->HasError()) {
			error = submitted->GetErrorObject();
		} else {
			stream = make_uniq<QueryResultStream<>>(std::move(submitted));
		}
	}

	bool HasError() const {
		return !stream || stream->HasError();
	}
	const string &GetError() const {
		return stream ? stream->GetError() : error.Message();
	}
	const ErrorData &GetErrorObject() const {
		return stream ? stream->GetErrorObject() : error;
	}
	//! Only meaningful when the query did not fail at submission
	ClientProperties GetClientProperties() const {
		return stream ? stream->GetClientProperties() : ClientProperties();
	}
	const vector<LogicalType> &GetTypes() const {
		return stream ? stream->GetTypes() : no_types;
	}
	const vector<Identifier> &GetNames() const {
		return stream ? stream->GetNames() : no_names;
	}
	unique_ptr<DataChunk> Fetch() {
		return stream ? stream->Fetch() : nullptr;
	}

private:
	unique_ptr<QueryResultStream<>> stream;
	ErrorData error;
	vector<LogicalType> no_types;
	vector<Identifier> no_names;
};

//! Submits a query and streams its result (v2.0 replacement for Connection::SendQuery)
inline unique_ptr<StreamingQuery> SubmitStream(Connection &conn, const string &sql) {
	return make_uniq<StreamingQuery>(conn.Submit(sql));
}

} // namespace duckdb
