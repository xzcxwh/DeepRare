# 操作系统相关功能：文件/目录操作
import os
# 文件操作：复制、移动、删除目录
import shutil
# 生成唯一任务 ID（通用唯一识别码）
import uuid
# 执行外部命令/脚本（运行变异分析流水线）
import subprocess
# 异步 IO 支持：FastAPI 后台任务
import asyncio
# 类型注解：字典、可选类型
from typing import Dict, Optional
# 日期时间：记录任务开始/结束时间
from datetime import datetime

# FastAPI Web 框架：构建 API 服务
from fastapi import FastAPI, UploadFile, File, BackgroundTasks, HTTPException
# 文件下载响应
from fastapi.responses import FileResponse, JSONResponse
# 数据模型验证：定义 API 返回格式
from pydantic import BaseModel

# ==================== 1. 初始化 FastAPI 应用 ====================
# 创建 API 实例，设置标题和版本
app = FastAPI(title="HERITA Variant Pipeline API", version="1.0.0")

# ==================== 2. 全局路径配置（项目目录结构） ====================
# 项目根目录：获取当前文件所在的顶层目录
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# 输入目录：存放用户上传的 VCF 文件
INPUT_DIR = os.path.join(PROJECT_ROOT, "input_data")
# 结果目录：流水线运行时的输出目录
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results")
# 归档目录：按任务 ID 保存历史结果（永久存储）
ARCHIVE_DIR = os.path.join(PROJECT_ROOT, "archived_results")
# 脚本目录：存放分析流水线脚本
SCRIPTS_DIR = os.path.join(PROJECT_ROOT, "scripts")

# ==================== 3. 自动创建必需目录 ====================
# 确保输入目录存在，不存在则创建
os.makedirs(INPUT_DIR, exist_ok=True)
# 确保归档目录存在
os.makedirs(ARCHIVE_DIR, exist_ok=True)


# ==================== 4. 任务状态定义 ====================
# 任务状态常量类
class TaskStatus:
    PENDING = "pending"  # 排队中：等待执行
    RUNNING = "running"  # 运行中：正在分析
    COMPLETED = "completed"  # 已完成：分析成功
    FAILED = "failed"  # 失败：任务出错


# API 任务信息数据模型（接口返回格式）
class TaskInfo(BaseModel):
    task_id: str  # 任务唯一 ID
    status: str  # 任务状态
    start_time: str  # 开始时间
    end_time: Optional[str] = None  # 结束时间（可选）
    message: str = ""  # 提示信息
    result_files: list = []  # 结果文件列表


# ==================== 5. 内存任务数据库 ====================
# 内存字典：存储所有任务的状态（服务重启后清空）
tasks: Dict[str, dict] = {}
# 异步锁：防止同时运行多个流水线（避免目录冲突）
pipeline_lock = asyncio.Lock()


# ==================== 6. 核心：运行分析流水线 ====================
def run_pipeline_task(task_id: str, input_vcf_path: str):
    """后台任务：执行基因变异分析流水线"""
    # 更新任务状态：运行中
    tasks[task_id]["status"] = TaskStatus.RUNNING
    tasks[task_id]["message"] = "Pipeline is running..."

    # 日志文件路径：按任务 ID 归档
    log_file = os.path.join(ARCHIVE_DIR, task_id, "pipeline.log")
    # 创建日志所在目录
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    try:
        # 构建执行命令：调用 run_pipeline.py 脚本
        cmd = [
            "python3",
            os.path.join(SCRIPTS_DIR, "run_pipeline.py"),
            input_vcf_path,  # 传入用户上传的 VCF 文件
            "--clean"  # 脚本参数：清理临时文件
        ]

        # 执行脚本，并将输出写入日志文件
        with open(log_file, "w") as f:
            process = subprocess.run(
                cmd,
                cwd=PROJECT_ROOT,  # 工作目录为项目根目录
                stdout=f,  # 标准输出写入日志
                stderr=subprocess.STDOUT,  # 错误信息也写入日志
                text=True  # 文本模式输出
            )

        # 命令执行成功（返回码=0）
        if process.returncode == 0:
            tasks[task_id]["status"] = TaskStatus.COMPLETED
            tasks[task_id]["message"] = "Pipeline completed successfully."

            # 归档结果：把运行结果复制到任务专属目录
            task_result_dir = os.path.join(ARCHIVE_DIR, task_id, "results")
            if os.path.exists(task_result_dir):
                shutil.rmtree(task_result_dir)
            shutil.copytree(RESULTS_DIR, task_result_dir)

            # 记录可下载的结果文件（VCF、报告、文本）
            tasks[task_id]["result_files"] = [
                f for f in os.listdir(task_result_dir)
                if f.endswith(".vcf") or f.endswith(".md") or f.endswith(".txt")
            ]
        else:
            # 执行失败
            tasks[task_id]["status"] = TaskStatus.FAILED
            tasks[task_id]["message"] = f"Pipeline failed with exit code {process.returncode}. Check logs."

    except Exception as e:
        # 代码异常报错
        tasks[task_id]["status"] = TaskStatus.FAILED
        tasks[task_id]["message"] = f"Internal error: {str(e)}"
    finally:
        # 无论成功/失败，记录结束时间
        tasks[task_id]["end_time"] = datetime.now().isoformat()


# ==================== 7. 异步包装器：带锁执行任务 ====================
async def process_pipeline(task_id: str, input_vcf_path: str):
    """包装器：获取锁并运行任务（防止并发冲突）"""
    async with pipeline_lock:  # 加锁：同一时间只运行一个任务
        loop = asyncio.get_event_loop()
        # 在线程池中运行同步任务，不阻塞 API
        await loop.run_in_executor(None, run_pipeline_task, task_id, input_vcf_path)


# ==================== 8. API 接口：启动流水线 ====================
@app.post("/pipeline/run", response_model=TaskInfo)
async def run_pipeline(
        background_tasks: BackgroundTasks,  # 后台任务
        file: UploadFile = File(...)  # 上传的 VCF 文件
):
    """上传 VCF 文件并启动变异分析流水线"""

    # 并发控制：已有任务在运行，直接返回错误
    if pipeline_lock.locked():
        raise HTTPException(status_code=429, detail="A pipeline task is already running. Please wait.")

    # 生成唯一任务 ID
    task_id = str(uuid.uuid4())
    # 时间戳：用于文件名
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # 保存上传文件到输入目录
    filename = f"{timestamp}_{task_id}_{file.filename}"
    input_path = os.path.join(INPUT_DIR, filename)

    # 写入文件
    with open(input_path, "wb") as buffer:
        shutil.copyfileobj(file.file, buffer)

    # 初始化任务状态
    tasks[task_id] = {
        "task_id": task_id,
        "status": TaskStatus.PENDING,
        "start_time": datetime.now().isoformat(),
        "message": "Task queued",
        "result_files": []
    }

    # 添加到后台任务（不阻塞前端）
    background_tasks.add_task(process_pipeline, task_id, input_path)

    # 返回任务信息
    return tasks[task_id]


# ==================== 9. API 接口：查询任务状态 ====================
@app.get("/pipeline/status/{task_id}", response_model=TaskInfo)
async def get_status(task_id: str):
    """查询任务运行状态（排队/运行/完成/失败）"""
    if task_id not in tasks:
        raise HTTPException(status_code=404, detail="Task not found")
    return tasks[task_id]


# ==================== 10. API 接口：获取运行日志 ====================
@app.get("/pipeline/logs/{task_id}")
async def get_logs(task_id: str):
    """获取任务运行日志（排查错误用）"""
    if task_id not in tasks:
        raise HTTPException(status_code=404, detail="Task not found")

    log_file = os.path.join(ARCHIVE_DIR, task_id, "pipeline.log")
    if not os.path.exists(log_file):
        return JSONResponse(content={"log": "Log file not created yet."})

    return FileResponse(log_file)


# ==================== 11. API 接口：下载结果文件 ====================
@app.get("/pipeline/download/{task_id}/{filename}")
async def download_result(task_id: str, filename: str):
    """下载分析结果文件（VCF / 报告文档）"""
    if task_id not in tasks:
        raise HTTPException(status_code=404, detail="Task not found")

    # 任务必须完成才能下载
    if tasks[task_id]["status"] != TaskStatus.COMPLETED:
        raise HTTPException(status_code=400, detail="Task not completed yet")

    file_path = os.path.join(ARCHIVE_DIR, task_id, "results", filename)
    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="File not found")

    return FileResponse(file_path, filename=filename)


# ==================== 12. 根路径接口 ====================
@app.get("/")
async def root():
    return {"message": "Welcome to HERITA Variant Pipeline API. Use /docs for API documentation."}