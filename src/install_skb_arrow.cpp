#include "install_skb_arrow.hpp"

#include "install_report.hpp"
#include "skb_arrow/install.hpp"

#include <atomic>
#include <string>

namespace duckdb {

namespace {

namespace sk = ::duckdb::miint::skb_arrow;

void InstallSkbArrowExecute(DataChunk &args, ExpressionState &state, Vector &result) {
	ClientContext &context = state.GetContext();
	sk::InstallOptions options;
	options.existing = sk::FindSkbArrow();
	options.uv = sk::FindUv();
	options.cache_dir = sk::MiintCacheDir();
	options.force = InstallForceArg(args);
	const sk::InstallReport report = sk::Install(options, [&context] { context.BlockingInterruptCheck(); });
	SetInstallReport(result, report.installed, report.path, report.version, report.message);
}

void SkbArrowAvailableExecute(DataChunk &args, ExpressionState &state, Vector &result) {
	// Cached once true, never while false: false is what a user sees before
	// running install_skb_arrow(), and it must not outlive the install. No
	// lock around the probe: a host that never answers holds up only the
	// query that started it, which can be interrupted.
	static std::atomic<bool> available {false};
	if (!available) {
		const std::string binary = sk::FindSkbArrow();
		ClientContext &context = state.GetContext();
		try {
			if (!binary.empty()) {
				sk::Probe(binary, [&context] { context.BlockingInterruptCheck(); });
				available = true;
			}
		} catch (const sk::Error &) {
			// Missing, broken, or incompatible: not available.
		}
	}
	result.Reference(Value::BOOLEAN(available), count_t(args.size()));
}

} // namespace

void RegisterSkbArrowInstallFunctions(ExtensionLoader &loader) {
	loader.RegisterFunction(InstallFunctionSet("install_skb_arrow", InstallSkbArrowExecute));
	// Volatile: a false answer folded into a plan would outlive the install.
	ScalarFunction available("skb_arrow_available", {}, LogicalType::BOOLEAN, SkbArrowAvailableExecute);
	available.SetVolatile();
	loader.RegisterFunction(available);
}

} // namespace duckdb
