#include "material/ScatteringMatrixReader.h"
#include <H5Cpp.h>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace phonomc {
namespace {

[[noreturn]] void invalid(const std::string& reason) {
    throw std::runtime_error("Scattering matrix HDF5: " + reason);
}

H5::DataSet dataset(const H5::H5File& file, const char* name) {
    if (H5Lexists(file.getId(), name, H5P_DEFAULT) <= 0)
        invalid(std::string("missing dataset ") + name);
    return file.openDataSet(name);
}

std::vector<hsize_t> dimensions(const H5::DataSet& data) {
    const auto space = data.getSpace();
    const int rank = space.getSimpleExtentNdims();
    if (rank < 0 || space.getSimpleExtentType() == H5S_NULL) invalid("invalid dataset dataspace");
    std::vector<hsize_t> shape(static_cast<std::size_t>(rank));
    if (rank > 0) space.getSimpleExtentDims(shape.data());
    return shape;
}

std::size_t check_shape(const H5::DataSet& data, const char* name, const std::vector<hsize_t>& expected) {
    if (dimensions(data) != expected) invalid(std::string(name) + " shape mismatch");
    std::size_t count = 1;
    for (const auto dim : expected) {
        if (dim > std::numeric_limits<std::size_t>::max() ||
            (count != 0 && dim > std::numeric_limits<std::size_t>::max() / count))
            invalid(std::string(name) + " dimensions overflow");
        count *= static_cast<std::size_t>(dim);
    }
    return count;
}

std::vector<unsigned long long> read_integers(const H5::DataSet& data, const char* name, std::size_t count) {
    if (data.getTypeClass() != H5T_INTEGER || data.getIntType().getPrecision() > 64)
        invalid(std::string(name) + " must contain integers of at most 64 bits");
    std::vector<unsigned long long> result(count);
    if (count == 0) return result;
    if (data.getIntType().getSign() == H5T_SGN_NONE) {
        data.read(result.data(), H5::PredType::NATIVE_ULLONG);
    } else {
        std::vector<long long> signed_values(count);
        data.read(signed_values.data(), H5::PredType::NATIVE_LLONG);
        for (std::size_t i = 0; i < count; ++i) {
            if (signed_values[i] < 0) invalid(std::string(name) + " must be nonnegative");
            result[i] = static_cast<unsigned long long>(signed_values[i]);
        }
    }
    return result;
}

std::vector<std::size_t> read_indices(const H5::H5File& file, const char* name, std::size_t count) {
    const auto data = dataset(file, name);
    check_shape(data, name, {static_cast<hsize_t>(count)});
    const auto raw = read_integers(data, name, count);
    std::vector<std::size_t> result(count);
    for (std::size_t i = 0; i < count; ++i) {
        if (raw[i] > std::numeric_limits<std::size_t>::max())
            invalid(std::string(name) + " index exceeds platform size range");
        result[i] = static_cast<std::size_t>(raw[i]);
    }
    return result;
}

std::vector<double> read_doubles(const H5::H5File& file, const char* name, const std::vector<hsize_t>& shape) {
    const auto data = dataset(file, name);
    const auto count = check_shape(data, name, shape);
    if (data.getTypeClass() != H5T_FLOAT && data.getTypeClass() != H5T_INTEGER)
        invalid(std::string(name) + " must be numeric");
    std::vector<double> result(count);
    if (count != 0) data.read(result.data(), H5::PredType::NATIVE_DOUBLE);
    return result;
}

} // namespace

ScatteringMatrixData read_scattering_matrix(const std::filesystem::path& path) {
    try {
        H5::H5File file(path.string(), H5F_ACC_RDONLY);
        const auto schema = dataset(file, "schema_version");
        check_shape(schema, "schema_version", {});
        if (read_integers(schema, "schema_version", 1)[0] != 1)
            invalid("unsupported schema_version (expected 1)");

        ScatteringMatrixData result;
        result.reference_temperature = read_doubles(file, "reference_temperature", {})[0];
        const auto modes = dataset(file, "mode_indices");
        const auto mode_shape = dimensions(modes);
        if (mode_shape.size() != 2 || mode_shape[1] != 2 || mode_shape[0] == 0 ||
            mode_shape[0] >= std::numeric_limits<std::size_t>::max() / 2)
            invalid("mode_indices must have nonempty shape (M, 2)");
        const auto n = static_cast<std::size_t>(mode_shape[0]);
        const auto mode_values = read_integers(modes, "mode_indices", n * 2);
        result.modes.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t axis = 0; axis < 2; ++axis) {
                const auto index = mode_values[2 * i + axis];
                if (index > static_cast<unsigned long long>(std::numeric_limits<int>::max()))
                    invalid("mode_indices exceeds integer mode range");
                result.modes[i][axis] = static_cast<int>(index);
            }
        }
        result.frequency_thz = read_doubles(file, "frequency", {mode_shape[0]});
        const auto values = dataset(file, "values");
        const auto value_shape = dimensions(values);
        if (value_shape.size() != 1 || value_shape[0] > std::numeric_limits<std::size_t>::max())
            invalid("values must have shape (nnz)");
        const auto nnz = static_cast<std::size_t>(value_shape[0]);
        result.column_offsets = read_indices(file, "column_offsets", n + 1);
        result.row_indices = read_indices(file, "row_indices", nnz);
        result.values = read_doubles(file, "values", value_shape);
        validate_scattering_matrix_data(result);
        return result;
    } catch (const H5::Exception& error) {
        invalid(path.string() + ": " + error.getDetailMsg());
    }
}

} // namespace phonomc
