const LOCAL_API_BASE_URL = "http://127.0.0.1:8000";

export function getApiBaseUrl(env = import.meta.env) {
  const configuredUrl = env?.VITE_API_BASE_URL?.trim();
  return (configuredUrl || LOCAL_API_BASE_URL).replace(/\/+$/, "");
}
