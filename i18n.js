const translations = {
  en: {
    nav_beg: "Beginner", nav_int: "Intermediate", nav_adv: "Advanced",
    welcome_title: "Welcome to the Stellar Structure and Evolution Simulator",
    welcome_desc: "The purpose of this simulator is for the user to be able to create and view simulations which evolve stars on a Hertzsprung-Russell (HR) diagram.",
    welcome_beg: "Archive of various stellar evolution simulations.",
    welcome_int: "Create simulations based on mass and metallicity.",
    welcome_adv: "Create simulations with more control parameters.",
    lbl_mass: "Mass of Stars (M☉)", lbl_z: "Metallicity Z", lbl_time: "Ending Time (Myr)",
    lbl_pts1: "pts1 (MS)", lbl_pts2: "pts2 (GB)", lbl_pts3: "pts3 (HG)", lbl_neta: "Reimers (neta)",
    btn_run: "Run Simulation", status_load: "Loading WASM engine…", status_ready: "Engine ready.",
    status_running: "Running simulation…", status_success: "Simulation successful!", status_error: "Error: ",
    col_tev: "Tev(Myr)", col_type: "Type", col_mt: "Mt", col_logl: "logL", col_logr: "logR", col_logt: "logT", col_spec: "Spec"
  },
  zh: {
    nav_beg: "初级", nav_int: "中级", nav_adv: "高级",
    welcome_title: "欢迎使用恒星结构与演化模拟器",
    welcome_desc: "本模拟器的目的是让用户能够创建并查看恒星在赫罗图 (HR Diagram) 上的演化模拟。",
    welcome_beg: "各种恒星演化模拟的存档。",
    welcome_int: "根据质量和金属丰度创建模拟。",
    welcome_adv: "使用更多控制参数创建模拟。",
    lbl_mass: "初始质量 (M☉)", lbl_z: "金属丰度 Z", lbl_time: "最大演化时间 (Myr)",
    lbl_pts1: "pts1 (MS)", lbl_pts2: "pts2 (GB)", lbl_pts3: "pts3 (HG)", lbl_neta: "Reimers 系数",
    btn_run: "开始模拟", status_load: "正在加载 WASM 引擎…", status_ready: "引擎已就绪。",
    status_running: "正在运行模拟…", status_success: "模拟成功！", status_error: "报错: ",
    col_tev: "Tev(Myr)", col_type: "类型", col_mt: "Mt", col_logl: "logL", col_logr: "logR", col_logt: "logT", col_spec: "光谱"
  }
};

let currentLang = new URLSearchParams(window.location.search).get('lang') || 'en';

function applyLang() {
  document.documentElement.lang = currentLang;
  document.querySelectorAll('[data-i18n]').forEach(el => {
    const key = el.getAttribute('data-i18n');
    if (translations[currentLang][key]) {
      el.textContent = translations[currentLang][key];
    }
  });
  const btn = document.getElementById('langBtn');
  if (btn) btn.textContent = currentLang === 'en' ? '中文' : 'English';
}

function toggleLang() {
  currentLang = currentLang === 'en' ? 'zh' : 'en';
  applyLang();
  // 同步更新 URL 方便分享
  const url = new URL(window.location);
  url.searchParams.set('lang', currentLang);
  window.history.pushState({}, '', url);
}

window.onload = applyLang;