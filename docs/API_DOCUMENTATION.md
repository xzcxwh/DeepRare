# HERITA Variant Pipeline API 文档

本文档描述了 HERITA 变异筛选流水线的 RESTful API 接口。该服务允许用户上传 VCF 文件，运行自动化分析流水线，并获取分析结果。

## 基本信息

- **服务名称**: HERITA Variant Pipeline API
- **版本**: 1.0.0
- **基础 URL**: `http://localhost:8000` (默认)
- **交互式文档**: `http://localhost:8000/docs` (Swagger UI)

## 接口列表

### 1. 提交分析任务

上传 VCF 文件并启动分析流水线。

- **URL**: `/pipeline/run`
- **方法**: `POST`
- **Content-Type**: `multipart/form-data`

#### 请求参数

| 参数名 | 类型 | 必选 | 描述 |
| :--- | :--- | :--- | :--- |
| `file` | File | 是 | 需要分析的 VCF 文件 (.vcf 或 .vcf.gz) |

#### 响应示例 (200 OK)

```json
{
  "task_id": "a829e080-f5b6-45a5-955b-922df36ff47d",
  "status": "pending",
  "start_time": "2025-12-08T19:02:21.317760",
  "end_time": null,
  "message": "Task queued",
  "result_files": []
}
```

#### 错误响应

- **429 Too Many Requests**: 当前已有任务在运行（服务仅支持单任务串行处理）。

---

### 2. 查询任务状态

获取指定任务的当前状态、进度信息和结果文件列表。

- **URL**: `/pipeline/status/{task_id}`
- **方法**: `GET`

#### 请求参数

| 参数名 | 位置 | 类型 | 必选 | 描述 |
| :--- | :--- | :--- | :--- | :--- |
| `task_id` | Path | string | 是 | 任务的唯一标识符 (UUID) |

#### 响应示例 (200 OK)

```json
{
  "task_id": "a829e080-f5b6-45a5-955b-922df36ff47d",
  "status": "completed",
  "start_time": "2025-12-08T19:02:21.317760",
  "end_time": "2025-12-08T19:02:43.446273",
  "message": "Pipeline completed successfully.",
  "result_files": [
    "final_intervar_classified.vcf",
    "STEP10_ACMG_DETAILED_LOG.txt",
    "STEP9_INTERVAR_ACMG_REPORT.md",
    "final_integrated.vcf"
  ]
}
```

#### 状态说明 (`status`)

- `pending`: 任务已创建，等待执行。
- `running`: 任务正在执行中。
- `completed`: 任务执行成功。
- `failed`: 任务执行失败。

---

### 3. 获取任务日志

获取任务执行过程中的详细日志输出（标准输出和标准错误）。

- **URL**: `/pipeline/logs/{task_id}`
- **方法**: `GET`

#### 请求参数

| 参数名 | 位置 | 类型 | 必选 | 描述 |
| :--- | :--- | :--- | :--- | :--- |
| `task_id` | Path | string | 是 | 任务的唯一标识符 |

#### 响应

- 返回纯文本格式的日志文件内容。

---

### 4. 下载结果文件

下载任务生成的特定结果文件。

- **URL**: `/pipeline/download/{task_id}/{filename}`
- **方法**: `GET`

#### 请求参数

| 参数名 | 位置 | 类型 | 必选 | 描述 |
| :--- | :--- | :--- | :--- | :--- |
| `task_id` | Path | string | 是 | 任务的唯一标识符 |
| `filename` | Path | string | 是 | 要下载的文件名（必须存在于 `result_files` 列表中） |

#### 常见结果文件

- `final_intervar_classified.vcf`: 最终包含 ACMG 分类的 VCF 文件。
- `final_integrated.vcf`: 整合了 ClinVar、dbscSNV 和 Top50 评分的 VCF 文件。
- `STEP10_ACMG_DETAILED_LOG.txt`: 详细的 ACMG 判定日志。
- `STEP9_INTERVAR_ACMG_REPORT.md`: ACMG 分类统计报告。

---

## 调用示例 (cURL)

### 提交任务

```bash
curl -X 'POST' \
  'http://localhost:8000/pipeline/run' \
  -H 'accept: application/json' \
  -H 'Content-Type: multipart/form-data' \
  -F 'file=@/path/to/your/input.vcf'
```

### 查询状态

```bash
curl -s http://localhost:8000/pipeline/status/YOUR_TASK_ID
```

### 下载最终结果

```bash
curl -O http://localhost:8000/pipeline/download/YOUR_TASK_ID/final_intervar_classified.vcf
```
