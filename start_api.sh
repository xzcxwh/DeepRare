#!/bin/bash
# 启动 HERITA API 服务

# 获取脚本所在目录
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_ROOT="$SCRIPT_DIR"

# 检查是否安装了依赖
if ! python3 -c "import fastapi" &> /dev/null; then
    echo "正在安装依赖..."
    pip3 install -r "$SCRIPT_DIR/api/requirements.txt"
fi

echo "启动 FastAPI 服务..."
echo "访问文档: http://localhost:8000/docs"

# 启动服务
# --reload 用于开发模式，生产环境可去掉
cd "$PROJECT_ROOT"
python3 -m uvicorn api.main:app --host 0.0.0.0 --port 8000 --reload
