// The COPY formats that plan through SequenceCopyPlan keep record order only when the query asks for it: with an ORDER
// BY, or with PRESERVE_ORDER true. Otherwise the records are a set, and the writer runs on all threads, as with
// PRESERVE_ORDER false. DuckDB's own default for a COPY is to keep insertion order.
//
// The SQL tests see the records a COPY writes, not the mode it ran in: a COPY that kept order when it need not is only
// slower, and one that ran on all threads can still write its records in order. These tests plan each COPY as DuckDB
// would run it, and read the mode DuckDB chose from the plan: PhysicalCopyToFile::parallel is true when the writer runs
// on all threads, and false when it keeps the order its rows come in.

#include <catch2/catch_test_macros.hpp>

#include "duckdb/execution/operator/persistent/physical_copy_to_file.hpp"
#include "duckdb/execution/physical_plan_generator.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/optimizer/optimizer.hpp"
#include "duckdb/parser/parser.hpp"
#include "duckdb/planner/planner.hpp"

#include <filesystem>
#include <string>

namespace {

using namespace duckdb;

// One output format: a table holding exactly the columns it writes, and the COPY options.
struct Output {
	std::string name;
	std::string table;
	std::string extension;
	std::string options;
};

const Output OUTPUTS[] = {
    {"FASTA", "fasta_rows", ".fasta", "FORMAT FASTA"},
    {"FASTQ", "fastq_rows", ".fastq", "FORMAT FASTQ"},
    {"SAM", "aln", ".sam", "FORMAT SAM, REFERENCE_LENGTHS 'ref_lens'"},
    {"headerless SAM", "aln", ".sam", "FORMAT SAM, INCLUDE_HEADER false"},
    {"BAM", "aln", ".bam", "FORMAT BAM, REFERENCE_LENGTHS 'ref_lens'"},
};

struct OrderFixture {
	DuckDB db;
	Connection con;

	OrderFixture() : db(nullptr), con(db) {
		for (auto sql :
		     {"SET threads = 8",
		      "CREATE TABLE fasta_rows AS SELECT 'r' || i AS read_id, 'ACGT' AS sequence1 FROM range(1000) t(i)",
		      "CREATE TABLE fastq_rows AS SELECT 'r' || i AS read_id, 'ACGT' AS sequence1, [30, 30, 30, "
		      "30]::UTINYINT[] AS qual1 FROM range(1000) t(i)",
		      "CREATE TABLE ref_lens AS SELECT 'ref' || k AS name, 1000000::BIGINT AS length FROM range(4) t(k)",
		      "CREATE TABLE aln AS SELECT 'r' || i AS read_id, 0::USMALLINT AS flags, 'ref' || (i % 4) AS reference, "
		      "(1 + i)::BIGINT AS position, 60::UTINYINT AS mapq, '4M' AS cigar, '*' AS mate_reference, 0::BIGINT AS "
		      "mate_position, 0::BIGINT AS template_length FROM range(1000) t(i)"}) {
			Run(sql);
		}
	}

	void Run(const std::string &sql) {
		auto result = con.Query(sql);
		INFO(sql << ": " << result->GetError());
		REQUIRE_FALSE(result->HasError());
	}

	// `COPY <source> TO <a path> (<output's options><extra>)`. Planning opens no file.
	static std::string Copy(const Output &output, const std::string &source, const std::string &extra = "") {
		auto target = std::filesystem::temp_directory_path() / ("miint_copy_order_never_written" + output.extension);
		return "COPY " + source + " TO '" + target.string() + "' (" + output.options + extra + ")";
	}

	// Plans `copy` as DuckDB would run it (bind, optimize, physical plan): true when its writer runs on all threads.
	bool RunsOnAllThreads(const std::string &copy) {
		INFO(copy);
		auto &context = *con.context;
		Parser parser(context);
		parser.ParseQuery(copy);
		REQUIRE(parser.statements.size() == 1);
		auto root_type = PhysicalOperatorType::INVALID;
		bool parallel = false;
		context.RunFunctionInTransaction([&]() {
			Planner planner(context);
			planner.CreatePlan(std::move(parser.statements[0]));
			Optimizer optimizer(*planner.binder, context);
			auto physical = PhysicalPlanGenerator(context).Plan(optimizer.Optimize(std::move(planner.plan)));
			auto &root = physical->Root();
			root_type = root.type;
			if (root_type == PhysicalOperatorType::COPY_TO_FILE) {
				parallel = root.Cast<PhysicalCopyToFile>().parallel;
			}
		});
		REQUIRE(root_type == PhysicalOperatorType::COPY_TO_FILE);
		return parallel;
	}
};

} // namespace

TEST_CASE("A sequence COPY with no ORDER BY runs on all threads", "[copy_order]") {
	OrderFixture fixture;
	for (auto &output : OUTPUTS) {
		INFO(output.name);
		// DuckDB keeps insertion order here by default (preserve_insertion_order is true), and the writers turn it off.
		CHECK(fixture.RunsOnAllThreads(OrderFixture::Copy(output, "(SELECT * FROM " + output.table + ")")));
		CHECK(fixture.RunsOnAllThreads(OrderFixture::Copy(output, output.table)));
	}
}

TEST_CASE("A sequence COPY whose query has an ORDER BY keeps its order", "[copy_order]") {
	OrderFixture fixture;
	for (auto &output : OUTPUTS) {
		INFO(output.name);
		const auto &table = output.table;
		CHECK_FALSE(fixture.RunsOnAllThreads(OrderFixture::Copy(output, "(FROM " + table + " ORDER BY read_id)")));
		// Under a projection, as a subquery or a CTE puts it.
		CHECK_FALSE(fixture.RunsOnAllThreads(
		    OrderFixture::Copy(output, "(SELECT * FROM (FROM " + table + " ORDER BY read_id) s)")));
		CHECK_FALSE(fixture.RunsOnAllThreads(
		    OrderFixture::Copy(output, "(WITH s AS (FROM " + table + " ORDER BY read_id) FROM s)")));
		// The optimizer makes ORDER BY with LIMIT a top-N.
		CHECK_FALSE(
		    fixture.RunsOnAllThreads(OrderFixture::Copy(output, "(FROM " + table + " ORDER BY read_id LIMIT 10)")));
	}
}

TEST_CASE("A sequence COPY keeps an ORDER BY when preserve_insertion_order is off", "[copy_order]") {
	OrderFixture fixture;
	fixture.Run("SET preserve_insertion_order = false");
	for (auto &output : OUTPUTS) {
		INFO(output.name);
		CHECK_FALSE(
		    fixture.RunsOnAllThreads(OrderFixture::Copy(output, "(FROM " + output.table + " ORDER BY read_id)")));
		CHECK(fixture.RunsOnAllThreads(OrderFixture::Copy(output, output.table)));
	}
}

TEST_CASE("PRESERVE_ORDER decides a sequence COPY's order, ORDER BY or not", "[copy_order]") {
	OrderFixture fixture;
	for (auto &output : OUTPUTS) {
		INFO(output.name);
		// true keeps the order the rows come in: a source's own order, as when converting a sorted file.
		CHECK_FALSE(fixture.RunsOnAllThreads(OrderFixture::Copy(output, output.table, ", PRESERVE_ORDER true")));
		// false drops an ORDER BY, as DuckDB documents.
		CHECK(fixture.RunsOnAllThreads(
		    OrderFixture::Copy(output, "(FROM " + output.table + " ORDER BY read_id)", ", PRESERVE_ORDER false")));
	}
}
