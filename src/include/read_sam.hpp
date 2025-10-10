#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "SAMRecord.hpp"

namespace duckdb {
class ReadSAMTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string sam_path;

		std::vector<std::string> names;            // field names
		std::vector<LogicalType> types;            // field types
		std::vector<miint::SAMRecordField> fields; // enum for convenience

		explicit Data(const std::string &sam_path)
		    : sam_path(sam_path),
		      names({"read_id",       "flags",           "reference", "position", "mapq",   "cigar",  "mate_reference",
		             "mate_position", "template_length", "tag_as",    "tag_xs",   "tag_ys", "tag_xn", "tag_xm",
		             "tag_xo",        "tag_xg",          "tag_nm",    "tag_yt",   "tag_md", "tag_sa"}),

		      types({LogicalType::VARCHAR,   // read_id
		             LogicalType::USMALLINT, // flags
		             LogicalType::VARCHAR,   // reference
		             LogicalType::INTEGER,   // position
		             LogicalType::UTINYINT,  // mapq
		             LogicalType::VARCHAR,   // cigar
		             LogicalType::VARCHAR,   // mate_reference
		             LogicalType::INTEGER,   // mate_positon
		             LogicalType::INTEGER,   // template length
		             LogicalType::INTEGER,   // tag_as
		             LogicalType::INTEGER,   // tag_xs
		             LogicalType::INTEGER,   // tag_ys
		             LogicalType::INTEGER,   // tag_xn
		             LogicalType::INTEGER,   // tag_xm
		             LogicalType::INTEGER,   // tag_xo
		             LogicalType::INTEGER,   // tag_xg
		             LogicalType::INTEGER,   // tag_nm
		             LogicalType::VARCHAR,   // tag_yt
		             LogicalType::VARCHAR,   // tag_md
		             LogicalType::VARCHAR}), // tag_sa

		      fields({miint::SAMRecordField::READ_ID,
		              miint::SAMRecordField::FLAGS,
		              miint::SAMRecordField::REFERENCE,
		              miint::SAMRecordField::POSITION,
		              miint::SAMRecordField::MAPQ,
		              miint::SAMRecordField::CIGAR,
		              miint::SAMRecordField::MATE_REFERENCE,
		              miint::SAMRecordField::MATE_POSITION,
		              miint::SAMRecordField::TEMPLATE_LENGTH,
		              miint::SAMRecordField::TAG_AS,
		              miint::SAMRecordField::TAG_XS,
		              miint::SAMRecordField::TAG_YS,
		              miint::SAMRecordField::TAG_XN,
		              miint::SAMRecordField::TAG_XM,
		              miint::SAMRecordField::TAG_XO,
		              miint::SAMRecordField::TAG_XG,
		              miint::SAMRecordField::TAG_NM,
		              miint::SAMRecordField::TAG_YT,
		              miint::SAMRecordField::TAG_MD,
		              miint::SAMRecordField::TAG_SA}) {
		}
	};
};

struct GlobalState : public GlobalTableFunctionState {
	mutex lock;
	miint::SAMReader reader;

	idx_t MaxThreads() const override {
		return 4;
	};

	GlobalState(const std::string &sam_path) : reader(miint::SAMReader(sam_path)) {};
};

static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
                                     vector<LogicalType> &return_types, vector<std::string> &names);

static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);

static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

static void SetResultVector(Vector &result_vector, const miint::SAMRecordField &field,
                            const std::vector<miint::SAMRecord> &records);
static void SetResultVectorNull(Vector &result_vector);
static void SetResultVectorString(Vector &result_vector, const miint::SAMRecordField &field,
                                  const std::vector<miint::SAMRecord> &records);
static void SetResultVectorStringNullable(Vector &result_vector, const miint::SAMRecordField &field,
                                          const std::vector<miint::SAMRecord> &records);
static void SetResultVectorListUInt8(Vector &result_vector, const miint::SAMRecordField &field,
                                     const std::vector<miint::SAMRecord> &records);

static TableFunction GetFunction();
};
}
; // namespace duckdb
