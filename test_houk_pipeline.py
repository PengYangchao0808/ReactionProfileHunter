import sys
import yaml
import logging
import glob
import shutil
import argparse
from pathlib import Path

# 添加项目根目录到路径
project_root = Path.cwd()
if str(project_root) not in sys.path:
    sys.path.append(str(project_root))

try:
    from rph_core.utils.qc_interface import run_gaussian_optimization, GaussianRouteBuilder, InputFactory
except ImportError as e:
    print(f"❌ 导入失败: {e}")
    sys.exit(1)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger("HoukTest")

def load_config(config_path):
    path = Path(config_path)
    if not path.exists():
        logger.error(f"❌ 找不到配置文件: {path}")
        sys.exit(1)
    with open(path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)

def main():
    parser = argparse.ArgumentParser(description="Houk Pipeline Final Verification")
    parser.add_argument("--config", default="config/defaults.yaml", help="Path to config file")
    args = parser.parse_args()

    print(">>> 🚀 Houk Pipeline Final Verification")
    print("-" * 60)

    # 1. 环境检查
    if not shutil.which("g16"):
        logger.error("❌ g16 命令未找到!")
        sys.exit(1)

    # 2. 加载配置
    config = load_config(args.config)
    logger.info(f"✅ 已加载配置: {args.config}")

    # 检查关键配置项
    if "dispersion" not in config["theory"]["optimization"]:
        logger.warning(">>> 没有提供 dispersion 配置，使用默认值 'GD3BJ'")
    if "solvent" not in config["theory"]["optimization"]:
        logger.warning(">>> 没有提供 solvent 配置，使用默认值 'acetone'")

    # 3. 验证 Route 构建
    print("\n>>> 🔄 调用核心模块构建 Route...")
    try:
        route = GaussianRouteBuilder.build(config)
        print(f">>> 📝 自动生成的 Route: [{route}]")
    except Exception as e:
        logger.error(f"❌ Route 构建失败: {e}")
        sys.exit(1)

    # 检查是否包含期望的关键字
    if "Def2SVP" in route and "em=GD3BJ" in route:
        print(">>> ✅ 验证通过: 翻译器工作正常 (Def2SVP & em=GD3BJ 均存在)")
    else:
        print(">>> ⚠️  警告: Route 可能缺少关键配置，请检查!")

    # 4. Dry-Run 预览
    res = config.get('resources', {})
    test_atoms = [
        {'symbol': 'O', 'x': 0.0, 'y': 0.0, 'z': 0.117},
        {'symbol': 'H', 'x': 0.0, 'y': 0.757, 'z': -0.469},
        {'symbol': 'H', 'x': 0.0, 'y': -0.757, 'z': -0.469}
    ]
    
    preview = InputFactory.create(route, test_atoms, 0, 1, res.get('mem','4GB'), res.get('nproc',4))
    print("\n>>> 📄 输入文件预览:")
    print("=" * 50)
    print(preview.strip())
    print("=" * 50)

    if input("\n>>> 确认执行? (y/n): ").lower() != 'y':
        sys.exit(0)

    # 5. 执行
    print("\n>>> 🔥 发射...")
    output_dir = Path("test_houk_final_output")
    output_dir.mkdir(exist_ok=True)
    
    result = run_gaussian_optimization(
        atoms=test_atoms,
        charge=0,
        mult=1,
        output_dir=output_dir,
        config=config,
        route=None 
    )

    if result['success']:
        print("\n>>> 🎉 成功 (SUCCESS)!")
        print(f">>> Log: {result['log_path']}")
    else:
        print(f"\n>>> ❌ 失败 (FAILURE): {result.get('error')}")

if __name__ == "__main__":
    main()
