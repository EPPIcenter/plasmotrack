from conan import ConanFile
from conan.tools.cmake import cmake_layout, CMakeToolchain

class ConanApplication(ConanFile):
    package_type = "application"
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeDeps"
    default_options = {
        "arrow/*:parquet": True,
        "arrow/*:with_zstd": True,
        "arrow/*:with_mimalloc": False,
        "thrift/*:with_openssl": False,
        "thrift/*:with_libevent": False,
    }

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.user_presets_path = False
        tc.generate()

    def requirements(self):
        # Linear algebra library
        self.requires("eigen/3.4.0")
        
        # Fast formatting library
        self.requires("fmt/9.1.0")
        
        # JSON library
        self.requires("nlohmann_json/3.11.3")
        
        # Boost C++ libraries (the arrow recipe requires >= 1.85)
        self.requires("boost/1.85.0")
        
        # Compression library
        self.requires("zlib/1.3.1")

        # Columnar storage: Arrow C++ 22 (the newest release that builds as C++17 with gcc >= 9) with Parquet and zstd
        self.requires("arrow/22.0.0")
