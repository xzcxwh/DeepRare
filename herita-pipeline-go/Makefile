.PHONY: build clean run test

# 默认目标
all: build

# 编译
build:
	@echo "🔨 编译 HERITA Pipeline..."
	@mkdir -p bin
	go build -o bin/pipeline ./cmd/pipeline
	@echo "✅ 编译完成: bin/pipeline"

# 清理
clean:
	@echo "🧹 清理..."
	rm -rf bin/
	@echo "✅ 清理完成"

# 运行（假设在 herita-project 根目录）
run: build
	@echo "🚀 运行流水线..."
	./bin/pipeline -root ..

# 运行测试
test:
	go test -v ./...

# 格式化代码
fmt:
	go fmt ./...

# 检查代码
vet:
	go vet ./...

# 下载依赖
deps:
	go mod download
	go mod tidy

# 交叉编译 Linux 版本
build-linux:
	@echo "🔨 编译 Linux 版本..."
	@mkdir -p bin
	GOOS=linux GOARCH=amd64 go build -o bin/pipeline-linux ./cmd/pipeline
	@echo "✅ 编译完成: bin/pipeline-linux"

# 帮助
help:
	@echo "HERITA Pipeline - Go 版本"
	@echo ""
	@echo "可用命令:"
	@echo "  make build       - 编译"
	@echo "  make run         - 编译并运行"
	@echo "  make clean       - 清理编译产物"
	@echo "  make test        - 运行测试"
	@echo "  make fmt         - 格式化代码"
	@echo "  make vet         - 检查代码"
	@echo "  make deps        - 下载依赖"
	@echo "  make build-linux - 交叉编译 Linux 版本"
