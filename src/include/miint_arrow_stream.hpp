#pragma once

#include "duckdb/common/arrow/arrow_appender.hpp"
#include "duckdb/common/arrow/arrow_converter.hpp"
#include "duckdb/common/arrow/arrow_wrapper.hpp"
#include "miint_streaming_query.hpp"

namespace duckdb {

//! Exposes a chunk source -- a StreamingQuery, or a materialized QueryResult -- as an ArrowArrayStream of
//! `batch_size`-row batches, converting on the consumer's thread.
//! Port of v1.5's ResultArrowArrayStreamWrapper, which DuckDB v2.0 removed (#25801). v2.0's replacement (ArrowFormat)
//! converts on worker threads and is not bounded by the streaming buffer, so miint keeps its own.
//! Ownership matches the old wrapper: the consumer's release() deletes this object.
template <class SOURCE>
class ChunkSourceArrowStream {
public:
	ArrowArrayStream stream;

	ChunkSourceArrowStream(unique_ptr<SOURCE> result_p, idx_t batch_size_p, ClientProperties client_properties_p)
	    : result(std::move(result_p)), batch_size(batch_size_p), client_properties(std::move(client_properties_p)) {
		if (batch_size == 0) {
			throw InvalidInputException("Arrow batch size must be greater than 0");
		}
		stream.private_data = this;
		stream.get_schema = GetSchema;
		stream.get_next = GetNext;
		stream.release = Release;
		stream.get_last_error = GetLastError;
	}

private:
	static ChunkSourceArrowStream &Self(ArrowArrayStream *s) {
		return *reinterpret_cast<ChunkSourceArrowStream *>(s->private_data);
	}

	static int GetSchema(ArrowArrayStream *s, ArrowSchema *out) {
		if (!s->release) {
			return -1;
		}
		auto &self = Self(s);
		out->release = nullptr;
		if (self.result->HasError()) {
			self.last_error = self.result->GetErrorObject();
			return -1;
		}
		try {
			ArrowConverter::ToArrowSchema(out, self.result->GetTypes(), IdentifiersToStrings(self.result->GetNames()),
			                              self.client_properties);
		} catch (std::exception &e) {
			self.last_error = ErrorData(e);
			return -1;
		}
		return 0;
	}

	static int GetNext(ArrowArrayStream *s, ArrowArray *out) {
		if (!s->release) {
			return -1;
		}
		auto &self = Self(s);
		out->release = nullptr;
		try {
			ArrowAppender appender(self.result->GetTypes(), self.batch_size, self.client_properties, {});
			idx_t count = 0;
			while (count < self.batch_size) {
				if (!self.current || self.offset >= self.current->size()) {
					self.current = self.result->Fetch();
					self.offset = 0;
					if (!self.current || self.current->size() == 0) {
						self.current.reset();
						break;
					}
				}
				auto take = MinValue(self.batch_size - count, self.current->size() - self.offset);
				appender.Append(*self.current, self.offset, self.offset + take, self.current->size());
				self.offset += take;
				count += take;
			}
			if (self.result->HasError()) {
				self.last_error = self.result->GetErrorObject();
				return -1;
			}
			if (count > 0) {
				*out = appender.Finalize();
			}
		} catch (std::exception &e) {
			self.last_error = ErrorData(e);
			return -1;
		}
		return 0;
	}

	static void Release(ArrowArrayStream *s) {
		if (!s || !s->release) {
			return;
		}
		s->release = nullptr;
		delete &Self(s);
	}

	static const char *GetLastError(ArrowArrayStream *s) {
		if (!s->release) {
			return "stream was released";
		}
		return Self(s).last_error.Message().c_str();
	}

	unique_ptr<SOURCE> result;
	idx_t batch_size;
	ClientProperties client_properties;
	unique_ptr<DataChunk> current;
	idx_t offset = 0;
	ErrorData last_error;
};

} // namespace duckdb
